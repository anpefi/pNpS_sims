#!/bin/bash

#SBATCH -J sel_03_McFL           # Job name
#SBATCH -o sel_03_McFL.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -a 1-100
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 02:00:00
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
NAME=sel_03_McFL_$SLURM_ARRAY_TASK_ID
$SCRIPT_DIR/run_simulation.R -i $NAME --seed 0 --nNS_pos 500 --nNS_neg 27000 \
          --nNS_neu 100 -c 1 --sampleEvery 0.0025 --reps 1 --model McFL \
          --finalTime 1000 --initSize 100000 --keepEvery 5 \
          --detectionSize 1e08
$SCRIPT_DIR/run_sampling.R -i $NAME --dir results


