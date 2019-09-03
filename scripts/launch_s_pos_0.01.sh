#!/bin/bash

#SBATCH -J s_pos_0.01           # Job name
#SBATCH -o s_pos_0.01.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 99:00:00
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
  --s_pos 0.01 --s_neg 0 \
  -r 100  -t 10000000 --keepEvery 10000000 \
  --n_pos 1 --n_neu 0 --n_neg 0 \
  --nNS_gene 500 --nS_gene 500 --prop_eff 1 \
  --initSize 20 --detectionSize 1e7 --mu 1e-8 \
  --wall.time 500000



