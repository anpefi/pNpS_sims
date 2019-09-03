#!/bin/bash

#SBATCH -J McF_50           # Job name
#SBATCH -o McF_50.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 5GB

set -e

## launch_simulations.sh

#load modules
module load gcc/6.4.0 R/2.6.0

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS.sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch
export RES_DIR=$PROJ_DIR/results

cd $RES_DIR
$SCRIPT_DIR/run_simulation2.R -i $SLURM_JOB_NAME --model McFL --sampleEvery 0.01 \
  --s_pos 0.25 --s_neg 0 \
  -r 100  -t 500000 --keepEvery 500000 \
  --n_pos 1 --n_neu 1 --n_neg 1 \
  --nNS_gene 500 --nS_gene 500 --prop_eff 0.5 \
  --initSize 20 --detectionSize 1e7 --mu 1e-8 \
  --wall.time 5000



