#!/bin/bash
#----------------------------------------------------
# Test for Uncertainty job script job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
# 
# Last revised: 4 feb 2021
#
# 
#----------------------------------------------------
#SBATCH -J create_cog_tiff           # Job name
#SBATCH -o cog.o 		# Name of stdout output file
#SBATCH -e cog.e       	# Name of stderr error file
#SBATCH -p normal	          	# Queue (partition) name
#SBATCH -N 1             				# Total # of nodes (must be 1 for serial)
#SBATCH -n 1            				# Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 5:30:00        				# Run time (hh:mm:ss)
#SBATCH --mail-user=wmobley@tamu.edu
#SBATCH --mail-type=all    				# Send email at begin and end of job
#SBATCH -A Flooding-and-Machine       	# Allocation name (req'd if you have more than 1)

eval "$(conda shell.bash hook)"
echo $CONDA_PREFIX

conda activate gdal
python /scratch/06659/wmobley/HDM/vrt.py

gdal_translate /scratch/06659/wmobley/HDM/outputFloodProb/damagePlain.vrt /scratch/06659/wmobley/HDM/outputFloodProb/damagePlain.tif -of COG -co COMPRESS=LZW -co BIGTIFF=YES -co NUM_THREADS=34 -co RESAMPLING=AVERAGE -co  OVERVIEW_RESAMPLING=AVERAGE

conda deactivate

