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


#SBATCH -J clip_hand           # Job name
#SBATCH -o hand.o       # Name of stdout output file
#SBATCH -e hand.e       # Name of stderr error file
#SBATCH -p skx-normal	          # Queue (partition) name

#SBATCH -N 3              # Total # of nodes (must be 1 for serial)
#SBATCH -n 15             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=wmobley@tamu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A PT2050-DataX       # Allocation name (req'd if you have more than 1)

module load launcher
pip install --user  pathlib urllib  contextlib
$LAUNCHER_DIR/paramrun

export LAUNCHER_JOB_FILE=hand_url.txt
export LAUNCHER_NJOBS=331 #number of jobs in the hand_url.txt
export LAUNCHER_PPN=5 #number of processes running on a node simultaneously
export LAUNCHER_NPROCS=15 #total number of processes running simultaneously
$LAUNCHER_DIR/paramrun

