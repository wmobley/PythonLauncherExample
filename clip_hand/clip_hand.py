import geopandas as gpd
from zipfile import ZipFile
import os
import pathlib
import sys
from contextlib import closing
import urllib2
import shutil
scratch = "/scratch/06659/wmobley/HDM/"

def unzip(zip_file):
    with ZipFile(os.path.join(scratch, zip_file), 'r') as zipObj:
       # Extract all the contents of zip file in current directory
       zipObj.extractall(os.path.join(scratch,zip_file.split(".")[0]))

def download_url(url, save_path):
    with closing(urllib2.urlopen(url)) as dl_file:
        with open(save_path, 'wb') as out_file:
            out_file.write(dl_file.read())

def get_tifs(Hand_url):

    zipped_file = Hand_url.split("/")[-1]
    huc = zipped_file.split(".")[0]
    if os.path.exists(os.path.join(scratch, "RawFiles/Hand", huc, "slope.tif" )):
        print(huc, "completed")
        return None
    download_url(Hand_url, os.path.join(scratch,zipped_file))
    print("downloaded: ", zipped_file)
    unzip(zipped_file)
    if os.path.exists(os.path.join(scratch,"RawFiles/Hand", huc))==False:
        os.makedirs(os.path.join(scratch, "RawFiles/Hand", huc))
    for in_file, out_file in [("", "elevation"), ("hand", "hand"), ("sd8", "slope")]:
        src = os.path.join(scratch,huc,  huc, huc+in_file+".tif")
        dst = os.path.join(scratch, "RawFiles/Hand", huc, out_file+".tif" )
        shutil.copyfile(src, dst, )
        print(dst)
    # shutil.rmtree(os.path.join(scratch,  huc))
    os.remove(os.path.join(scratch, zipped_file))
    print("finished: ", huc)
for arg in sys.argv[1:]:
    get_tifs(arg)