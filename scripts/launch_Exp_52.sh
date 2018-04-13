#!/bin/bash

#SBATCH -J Exp_52           # Job name
#SBATCH -o Exp_52.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -a 1-400
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu 20GB

set -e

## launch_simulations.sh

#load modules
module load gcc R

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS.sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch
export RES_DIR=$PROJ_DIR/results

cd $RES_DIR
NAME=sel_52_Exp_$SLURM_ARRAY_TASK_ID
$SCRIPT_DIR/run_simulation2.R -i $NAME --s_pos 1 --s_neg -0.1 -r 25

