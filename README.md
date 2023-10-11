# PythonLauncherExample
Examples of Sbatch scripts I've commonly used for Python for Lonestar 6. 

This Repository has three examples, One that was used to clip datasets by watershed; one that clipped hand datasets to smaller watersheds, and one that creates a large Cloud Optimized Geotif. 

## clip_Dataset
THis directory has three files in it. The Sbatch sets up the parallelized python script for clip_datasets.py through the clip_datset_list.txt file. 
This provides and example of how to use the launcher module on 1 node with multiple processes running simultaneously. 

## clip_hand
This directory shows how to use launcher across multiple nodes. 

## CreateCOG
This directory provides code to generate a large cloud optimized geotif. 

