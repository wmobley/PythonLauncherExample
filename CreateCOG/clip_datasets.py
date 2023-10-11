import os
import pandas as pd
import geopandas as gpd
import gc
import numpy as np
import rasterio
import rasterio.mask
import shutil
import os
from Hazard_Estimates.Model import *
import argparse
from rasterio.warp import calculate_default_transform, reproject, Resampling
from scipy import ndimage
import sys
import gdal
import richdem as rd
from multiprocessing import cpu_count, Pool
import time




roughness = {21: 0.0404,

                22: 0.0678,
                23: 0.0678,
                24: 0.0404,
                31: 0.0113,
                41: 0.36,
                42: 0.32,
                43: 0.40,
                52: 0.40,
                71: 0.368,
                81: 0.325,
                90: 0.086,
                95: 0.1825}

def rescale(root, name):
        # This Function rescales the given rasters based on the highest resolution image used. Currently that is the hand.tif at 10m.
        with rasterio.open(os.path.join(root, 'hand.tif')) as mask:
            # with rasterio.open(os.path.join(root, 'demfill.tif')) as mask2:
            shape = [{'type': 'Polygon', 'coordinates': [[(mask.bounds.left, mask.bounds.top),
                                                        (mask.bounds.left, mask.bounds.bottom),
                                                        (mask.bounds.right, mask.bounds.bottom),
                                                        (mask.bounds.right, mask.bounds.top)]]}]

            with rasterio.open(os.path.join(root, "temp.tif")) as src:  # loads a tempory tif previously clipped.
                transform, width, height = calculate_default_transform(
                    # note temp.tif is developed in the clip to watershed function below.
                    src.crs, mask.crs, mask.width, mask.height, *mask.bounds)
                kwargs = src.meta.copy()
                kwargs.update({
                    'crs': mask.crs,
                    'transform': transform,
                    'width': width,
                    'height': height,

                })
                with rasterio.open(os.path.join(root, name), 'w',
                                **kwargs) as dst:  # save the temporary file now rescaled to 10-meters
                    for i in range(1, src.count + 1):
                        reproject(
                            source=rasterio.band(src, i),
                            destination=rasterio.band(dst, i),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=mask.crs,
                            resampling=Resampling.bilinear)
        os.remove(os.path.join(root, "temp.tif"))


def safely_reduce_dtype(ser):  # pandas.Series or numpy.array
    # reduces the datatype to the lowest possible to reduce storage.
    orig_dtype = "".join([x for x in ser.dtype.name if x.isalpha()])  # float/int
    mx = 1

    new_itemsize = np.min_scalar_type(ser).itemsize
    if mx < new_itemsize:
        mx = new_itemsize
    new_dtype = orig_dtype + str(mx * 8)
    return new_dtype


def clip_to_boundary(in_directory, out_directory, boundary_geom,
                    in_raster, out_raster):
    # Clips the in raster based on the boundary geometry. Usually a Hydrologic Unit.
    if os.path.exists(os.path.join(out_directory, out_raster)): return None
    with rasterio.open(os.path.join(in_directory, in_raster)) as src:
        try:
            out_image, out_transform = rasterio.mask.mask(src, boundary_geom, crop=True)
        except:
            print("error", out_directory, out_raster)
        out_meta = src.meta
        
        if out_raster in [f"Distance2Lakes.tif", f"Distance2Streams.tif",f"Distance2Coast.tif"]:
            nodata = 0
        else:
            nodata = src.nodata
        out_meta.update({"driver": "GTiff",
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "transform": out_transform,
                        "dtype": out_image.dtype,
                        "nodata":nodata}
                        )
        if out_raster == "hand.tif":  # if the raster is a hand raster then it is the highest resolution so write the image otherwise rescale the image.
            with rasterio.open(os.path.join(out_directory, out_raster), 'w', **out_meta) as dest:
                dest.write(out_image)
        else:
            with rasterio.open(os.path.join(out_directory, "temp.tif"), 'w', **out_meta) as dest:
                dest.write(out_image)
            rescale(out_directory, out_raster)


def clip_roughness(directory, boundary, year):
    # calcuate roughness coefficient based on landuse
    clip_to_boundary("RawFiles/Landcover", directory, boundary,  # clip the landcover to the boundary
                    f"NLCD_{year}_Land_Cover_L48_20190424.img",
                    f"Landcover{year}.tif")
    with rasterio.open(os.path.join(directory,
                                    f"Landcover{year}.tif")) as src:  # load landcover and use a lookup table to estimate roughness.
        image = src.read(1)
        u, inv = np.unique(image, return_inverse=True)
        img = np.array([roughness.get(x, 0) for x in u])[inv].reshape(image.shape)

        out_meta = src.meta

        out_meta.update({"driver": "GTiff",
                        "height": img.shape[0],
                        "width": img.shape[1],
                        "dtype": "float32",
                        })

        with rasterio.open(os.path.join(directory, f"roughness{year}.tif"), 'w', **out_meta) as dst:
            dst.write(img.astype("float32"), 1)
    os.remove(os.path.join(directory, f"Landcover{year}.tif"))


def clip_twi(directory):
    # Calculate TWI based on slope and flow accumulation.
    TWI_Save = "{}/TWI.tif".format(directory)
    if os.path.exists(TWI_Save): return None
    with rasterio.open(os.path.join(directory, f"FlowAccumulation.tif")) as src_acc:
        flow_accumulation = src_acc.read(1)
        with rasterio.open(os.path.join(directory, f"slope.tif")) as src_slope:
            slope = src_slope.read(1)
            img = np.log(((flow_accumulation * 900) + 1) / (np.tan((slope + 0.000001) / (180 / np.pi))))
            out_meta = src_acc.meta

            out_meta.update({"driver": "GTiff",
                            "height": img.shape[0],
                            "width": img.shape[1],
                            "dtype": img.dtype,
                            })
            with rasterio.open(TWI_Save, 'w', **out_meta) as dst:
                dst.write(img, 1)
def weighted_accum(in_dir, weight, out_directory, out_raster, boundary_geom, dem):
            # uses the pysheds library to calculate flow accumulation, however its weighted with by a given array.
            if os.path.exists(os.path.join(out_directory, out_raster)): return None
            with rasterio.open(os.path.join(out_directory, weight)) as src:
                weights = src.read(1)
                weights = np.where(weights == -9999, 0, 1)
                norm = np.linalg.norm(weights)
                weights = weights / norm
                #
                # g.accumulation(data='dir', weights=weights, out_name='weights_accum')
                # grid.to_raster('weights_accum', os.path.join(out_directory, out_raster),dtype=np.int32)

                accum = rd.FlowAccumulation(dem, method='D8', weights=weights.astype('float64'))

                rd.SaveGDAL(os.path.join(out_directory, "temp2.tif"), accum)
                print(np.max(accum))
                clip_to_boundary(out_directory, out_directory, boundary_geom, f"temp2.tif",
                                out_raster)
                os.remove(os.path.join(out_directory, "temp2.tif"))




def clip(huc12):

        os.chdir("/mnt/corral-sync/HDM/web.corral.tacc.utexas.edu/HDM/")
        
    # Loop through each polygon in shapefile.
        
        shapefile = gpd.read_file("/mnt/code-391ff5ac-6576-460f-ba4d-7e03433c68b6/Users/wmobley/Code/images/huc12_damageplain.shp")
        # shapefile = shapefile.loc[(shapefile.HUC12==huc12)]  # running one of these currently
        shapefile.reset_index(inplace=True)

        geom = shapefile.geometry.values
        huc12 = shapefile.huc12.values[0]
        
        print(huc12)
       
        
        dst_crs = 'EPSG:5070'
        
        
        out_dir = f"Spatial_Index/huc{huc12[:-4]}/huc{huc12}"
       
       
        
        size=[]
        if os.path.exists(f"Spatial_Index/huc{huc12[:-4]}") == False:
            try:
                 os.makedirs(f"Spatial_Index/huc{huc12[:-4]}")
            except FileExistsError as e:
                print(e)
        if os.path.exists(out_dir) == False:
            print(out_dir)
            os.makedirs(out_dir)  # check to see if this directory exists if not make it.
    
        
        for topography in ['hand', 'elevation']:
            if os.path.exists(os.path.join(f"RawFiles/Hand/{huc12[:-6]}", f"{topography}_proj.tif")) == False:
                print(os.path.join(f"RawFiles/Hand/{huc12[:-6]}", f"{topography}_proj.tif"))
                in_hand = os.path.join(f"RawFiles/Hand/{huc12[:-6]}", f"{topography}.tif")
                out_hand = os.path.join(f"RawFiles/Hand/{huc12[:-6]}", f"{topography}_proj.tif")
                with rasterio.open(in_hand) as src:
                    transform, width, height = calculate_default_transform(
                        src.crs, dst_crs, src.width, src.height, *src.bounds)
                    kwargs = src.meta.copy()
                    kwargs.update({
                        'crs': dst_crs,
                        'transform': transform,
                        'width': width,
                        'height': height
                    })

                    with rasterio.open(out_hand, 'w', **kwargs) as dst:
                        for i in range(1, src.count + 1):
                            reproject(
                                source=rasterio.band(src, i),
                                destination=rasterio.band(dst, i),
                                src_transform=src.transform,
                                src_crs=src.crs,
                                dst_transform=transform,
                                dst_crs=dst_crs,
                                resampling=Resampling.bilinear)

            if topography == "elevation":
                clip_to_boundary(
                    f"RawFiles/Hand/{huc12[:-6]}",
                    out_dir, geom, f"{topography}_proj.tif", f"dem.tif")  # Clip Hand
            else:
                clip_to_boundary(
                    f"RawFiles/Hand/{huc12[:-6]}",
                    out_dir, geom, f"{topography}_proj.tif", f"{topography}.tif")  # Clip Hand
        #This for loop calculates euclidean distance.
        for water in ['Lakes', "Coast", "Streams"]:
                    clip_to_boundary("RawFiles", out_dir, geom.buffer(200), f"Distance2{water}_.tif",
                            f"Distance2{water}.tif")
        
        clip_to_boundary("RawFiles/Topography", out_dir, geom, f"elevation.tif",
                         f"dem30m.tif")

        # clip_to_boundary("RawFiles/Topography", out_dir, geom, f"texas_slope.tif",
        #                  f"slope30m.tif")
        gc.collect()  # clean up ram
        in_elevation = os.path.join(out_dir, f"dem.tif")
        dem = rd.LoadGDAL(in_elevation)
        rd.FillDepressions(dem, epsilon=True, in_place=True)
        slope = rd.TerrainAttribute(dem, attrib='slope_riserun')
        rd.SaveGDAL(os.path.join(out_dir, 'slope.tif'), slope)
        time.sleep(0.5)
        in_elevation = os.path.join(out_dir, f"dem30m.tif")
        dem = rd.LoadGDAL(in_elevation)
        accum_d8 = rd.FlowAccumulation(dem, method='D8')
        # rd.SaveGDAL(os.path.join(out_dir, 'FlowAccumulation.tif'), accum_d8)

        # Once slope and flow acculation are clipped then TWI can be calculated.
        clip_twi(out_dir)
        # This clips rainfall intensities for specific storms. Not necessary for the first analysis but needed later down the line.
        for hr in [1, 2, 3, 4, 8, 12, 24, ]:
            for storm in [
                'taxday',
                'harvey']:
                clip_to_boundary(r"RawFiles/Precip/{}/intensity/projected".format(storm), out_dir, geom,
                                 f"{storm}{hr}hr.tif",
                                 f"{storm}{hr}hr.tif")

        # Clip KSAT and then generate the accumulated KSAT.

        
        clip_to_boundary(r"RawFiles", out_dir, geom, "Ksat.tif",
                        f"ksat.tif")
        weighted_accum(out_dir,
                    "ksat.tif",
                    out_dir, 'AverageKSAT.tif', geom, dem)
            
        # These probabilistic precipitations. 12hr and 60 minutes. We use three probabilities 25, 100, and 500 year.
        for year in [25, 100, 500]:
            for hour in ["12ha", "60ma"]:
                print(out_dir,f"Precip{year}yr_{hour}_cog.tif")
                try:
                    clip_to_boundary("RawFiles/Precip", out_dir, geom, f"Precip{year}yr_{hour}_cog.tif",
                                f"Precip{year}_{hour}.tif")
                except ValueError as e:
                    print(huc, e)
                    with open('/home/azureuser/cloudfiles/code/Users/wmobley/Code/hucErrors.txt', 'a') as the_file:
                            the_file.write(f"{huc}, Precip{year}_{hour}.tif, \n" )
        #Thesea are the dynamic rasters. Using Imperviousness and Landcover.
        #It iterates over each year, and then clips a given raster. Impervious 2016 was named different so it stands along
        for i in [
            2001, 2004, 2006,
            2008, 2011, 2013,
            2016]:
            if i == 2016:
                clip_to_boundary("RawFiles/Impervious", out_dir, geom, f"impervious2016.tif",
                                f"impervious{i}.tif")
            elif i in [2001, 2006, 2011]:
                clip_to_boundary(r"RawFiles/Impervious/nlcd_{}_impervious_2011_edition_2014_10_10".format(i),
                                out_dir, geom,
                                f"nlcd_{i}_impervious_2011_edition_2014_10_10.img",
                                f"impervious{i}.tif")
         
            if i in [2001,2006,2011,2016]:
                time.sleep(0.5)
                try:
                    clip_roughness(out_dir, geom, i)
                except:
                    print(f"roughnessError {huc} {i}")
                    with open('hucErrors.txt', 'a') as the_file:
                            the_file.write(f"roughnessError {huc}" )
                time.sleep(0.5)
                try: 
                    weighted_accum(out_dir,
                            f"roughness{i}.tif",
                            out_dir, f'AverageRoughness{i}.tif', geom, dem)
                except:
                    print(f"average Roughness error {huc}")
                    with open('hucErrors.txt', 'a') as the_file:
                            the_file.write(f"average Roughness error {huc} {i}" )
   
def make_VRT(args):

    file_ =args.get("file")
    directory =  args.get("directory")
    sub_directories = [os.path.join(directory, name, file_) for name in os.listdir(directory) if
                    os.path.isdir(os.path.join(directory, name))]
    
    with rasterio.open(sub_directories[0])as src:
        srcnodata = src.nodatavals[0]
        if srcnodata==None:
            srcnodata=0
       
        print(f"{file_[:-4]}.vrt")
        gdal.BuildVRT(os.path.join(directory, f"{file_[:-4]}.vrt"), sub_directories,options=gdal.BuildVRTOptions(srcNodata=srcnodata, VRTNodata=srcnodata) )



if __name__=="__main__":
   
    
       
    sfile = gpd.read_file(r"./images/huc12_damageplain.shp")
    

    if len(sys.argv)>1:
        
        hucs = sys.argv[1:]
    else:
        hucs = sfile.huc12.apply(lambda row: row).unique()
    pool = Pool(cpu_count()-1)
    # for huc in hucs:
    #     # huc12 = sfile.loc[(sfile.HUC12.str.startswith(huc))].HUC12.values  # running one of these currently
        
        
    #     
    #     # pool = Pool(4)
    #     print(cpu_count())
    #     # pool.map(clip, huc12)
    #     [clip(h) for h in sfile.huc12]
    files= []
#  sfile.HUC12.apply(lambda row: row[:8]).unique()
    for huc in hucs:
        directory = f"/mnt/corral-sync/HDM/web.corral.tacc.utexas.edu/HDM/Spatial_Index/huc{huc[:8]}"
        
        # for year in[2001,2011,2006,2016]:
            
            
        #    files.append({"file":f"roughness{year}.tif","directory":directory})
            
        #    files.append({"file":f'impervious{year}.tif',"directory":directory})
        #    files.append({"file":f'AverageRoughness{year}.tif',"directory":directory})
        # for hr in [1, 2, 3, 4, 8, 12, 24, ]:
        #     for storm in [
        #         'taxday',
        #         'harvey']:
        #         files.append({"file":f"{storm}{hr}hr.tif","directory":directory})
               
        # for year in [25, 100, 500]:
        #    for hour in ["12ha", "60ma"]:
        #        files.append({"file":f"Precip{year}_{hour}.tif","directory":directory})
        # for water in ['Lakes', "Coast", "Streams"]:
        #     files.append({"file":f"Distance2{water}.tif","directory":directory})
        for listed_file in [
            # f"FlowAccumulation.tif",
            # f"slope.tif",
        #    f"dem.tif",
           'AverageKSAT.tif',
        #    'hand.tif'
        ]:
           files.append({"file":listed_file,"directory":directory})
    
    # pool.map( make_VRT, files)
    make_VRT({"file":f"AverageKSAT.tif","directory":"/mnt/corral-sync/HDM/web.corral.tacc.utexas.edu/HDM/  t tSpatial_Index"})