#!/bin/bash

#SBATCH -J slim_pos_0.5           # Job name
#SBATCH -o slim_pos_0.5.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu 15GB  ## Set 16 for cancer

set -e

## launch_simulations.sh

#load modules
module load slim

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS.sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch
export RES_DIR=$PROJ_DIR/results

cd $RES_DIR
for r in {1..100}
do
    slim -d s=0.5 -d mu=1e-8 -d Nmax=1e7 -d tmax=1000000 -d "name='$SLURM_JOB_NAME'" $SCRIPT_DIR/McFarland.slim 
done


