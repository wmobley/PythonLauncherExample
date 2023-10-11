#!/bin/bash
#----------------------------------------------------
# Test for Uncertainty job script job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
# 
# Last revised: 9 sept 2019
#
# 
#----------------------------------------------------


#SBATCH -J clip_watersheds           # Job name
#SBATCH -o watersheds.o # Name of stdout output file
#SBATCH -e watersheds.e       # Name of stderr error file
#SBATCH -p skx-normal	          # Queue (partition) name

#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 04:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=wmobley@tamu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Flooding-and-Machine       # Allocation name (req'd if you have more than 1)

eval "$(conda shell.bash hook)"
echo $CONDA_PREFIX

conda activate gdal
pip install --upgrade git+https://github.com/wmobley/Hazard_Estimates.git
module load launcher

$LAUNCHER_DIR/paramrun

export LAUNCHER_JOB_FILE=clip_dataset_list.txt
export LAUNCHER_NJOBS=45 #Total Number of jobs in your job file
export LAUNCHER_PPN=45 #total number of processes you want per node
export LAUNCHER_NPROCS=45 #Total Number of Processes you want running simulataneously. 
$LAUNCHER_DIR/paramrun
conda deactivate

