#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 30GB
#SBATCH -t 02:00:00

#load modules
module load gcc R

## sampling.sh
set -e
set -u

## PARSING ARGUMENTS
while getopts 'i:d:h' OPTION; do
  case "$OPTION" in
    i)
      ID="$OPTARG"
      ;;

    d)
      DL="$OPTARG"
      ;;

    h)
      echo "Help don't available"
      ;;
    ?)
      echo "script usage: $(basename $0) [-i ID] [-d detectio_limit] [-h]" >&2
      exit 1
      ;;
  esac
done

# DEFINE ENVIRONMENT VARIABLES
export PROJ_DIR=$HOME/pNpS_sims
export SCRIPT_DIR=$PROJ_DIR/scripts
export SCRATCH=$PROJ_DIR/scratch

cd $SCRATCH
$SCRIPT_DIR/run_sampling.R -i $ID -d $DL