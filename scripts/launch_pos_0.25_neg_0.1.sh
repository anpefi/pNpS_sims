#!/bin/bash

#SBATCH -J pos_0.25_neg_0.1           # Job name
#SBATCH -o out/pos_0.25_neg_0.1_%j    # Specify stdout output file (%j expands to jobId)
#SBATCH -n 1                     # Total number of tasks
#SBATCH -t 96:00:00
#SBATCH --mem-per-cpu 10GB

set -e

## launch_simulations.sh

#load modules
module load gcc/6.4.0 R/3.6.0

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS.sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch
export RES_DIR=$PROJ_DIR/results

cd $RES_DIR
$SCRIPT_DIR/run_simulation3.R -i ${SLURM_JOB_NAME}_${SLURM_JOB_ID} --model Exp --sampleEvery 0.01 \
  --s_pos 0.25 --s_neg -0.1 \
  -r 1000  -t 10000000 --keepEvery 10000000 \
  --n_pos 1 --n_neu 1 --n_neg 1.1 \
  --nNS_gene 5000 --nS_gene 5000 --prop_eff 0.002 \
  --initSize 20 --detectionSize 1e4 --mu 1e-6 \
  --wall.time 50000



