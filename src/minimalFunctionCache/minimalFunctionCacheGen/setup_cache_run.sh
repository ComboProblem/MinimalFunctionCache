#!/bin/bash -x

echo "Running on $(hostname)"
echo "Setting up run."
source ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh

module load apptainer
apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/cache_job_setup.py $NUM_BKPT $SAMPLE_SIZE $TIME_PER_BATCH $MAX_NUM_ROW $MAX_STD $BACKEND

ls ~/MinimalFunctionCache/TEMP
