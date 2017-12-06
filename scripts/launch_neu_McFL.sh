#!/bin/bash

#SBATCH -J neu_McFL           # Job name
#SBATCH -o neu_McFL-%a.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -a 1-100
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 01:00:00
#SBATCH --mem-per-cpu 6GB


## launch_simulations.sh

#load modules
module load gcc R

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS_sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch
export RES_DIR=$PROJ_DIR/results

cd $RES_DIR
NAME=neu_McFL_$SLURM_ARRAY_TASK_ID
$SCRIPT_DIR/run_simulation.R -i $NAME --seed 0 --nNS_pos 0 --nNS_neg 0 \
          --nNS_neu 27600 -c 1 --sampleEvery 0.01 --reps 10 --model McFL \
          --finalTime 50000 --keepEvery 500 --detectionSize 30000
$SCRIPT_DIR/run_sampling.R -i $NAME -d 0.000000001
$SCRIPT_DIR/run_sampling.R -i $NAME -d 0.001
$SCRIPT_DIR/run_sampling.R -i $NAME -d 0.01
$SCRIPT_DIR/run_sampling.R -i $NAME -d 0.05
