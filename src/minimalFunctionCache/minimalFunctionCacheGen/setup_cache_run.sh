#!/bin/bash -x
export K="$1"

echo "Running on $(hostname)"
echo "Setting up run for $K breakpoints."

source ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
module load apptainer
apptainer run $MFC_TARGET/cgf.sif python3 ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_job_setup.py 
