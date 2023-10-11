import os
import gdal
import fiona
import rasterio
from rasterio.mask import mask
gdal.SetConfigOption('GDAL_MAX_DATASET_POOL_SIZE', '7000')
def make_VRT(directory, files):
    gdal.BuildVRTOptions(VRTNodata="nan")
    gdal.BuildVRT(os.path.join(directory, f"damagePlain.vrt"), files)

directory = "/mnt/corral-sync/HDM/web.corral.tacc.utexas.edu/HDM/outputFloodProb"
for root, dirs, files in os.walk(directory, topdown=True):
   
    # if os.path.exists(os.path.join(root, f"damagePlain.tif")): continue
    vrt_array = []
    if(len(files))==0:
        continue
    for f in files:
        
        if f.endswith(".tif"):
           
            if  f.endswith("MultiModel.tif") or f.endswith("noweight.tif") or f == "damagePlain.tif" or  f == "damagePlain_State.tif" or root.endswith("clipped"):continue
   
            vrt_array.append(os.path.join(root,f))
   
    make_VRT(root, vrt_array)

    gdal.TranslateOptions( format ="COG",
                           )
    gdal.Translate(os.path.join(root, f"damagePlain.tif"),
                   os.path.join(root, f"damagePlain.vrt")
                   )
    print(root)
print("HUC damage plains")
vrt_array = []
for root, dirs, files in os.walk(directory, topdown=True):


    if(len(files))==0:
        continue
    for file in files:
        if file=="damagePlain.tif":
            vrt_array.append(os.path.join(root,file))

make_VRT(directory, vrt_array)
gdal.Translate(os.path.join(directory, f"damagePlain.tif"),
                   os.path.join(directory, f"damagePlain.vrt")
                   )
