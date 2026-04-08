#!/bin/bash -x
export K="$1"
export BKPTS_PATH="$2"
export REP_ELEMS_PATH="$3"
echo "Running on $(hostname)"
echo "Setting up run."
echo "Working on $K"

source ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
module load apptainer
apptainer run cgf.sif python3 ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_job_setup.py 
