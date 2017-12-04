#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 30GB
#SBATCH -t 02:00:00
## launch_simulations.sh

#load modules
module load gcc R

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS_sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch


cd $SCRATCH
$SCRIPT_DIR/run_simulation.R -i neu_Exp --seed 1234 --nNS_pos 0 --nNS_neg 0 --nNS_neu 27600 -c 4 --reps 100 --model Exp
