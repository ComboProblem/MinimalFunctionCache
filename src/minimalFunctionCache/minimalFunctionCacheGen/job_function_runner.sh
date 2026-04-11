#!/bin/bash

export K="$1"

echo "Running on $(hostname)"
echo "Task: $SLURM_ARRAY_TASK_ID"
source ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
source ~/MinimalFunctionCache/TEMP/temp_job_info_$K.sh
module load apptainer
apptainer run $MFC_TARGET/cgf.sif python3 ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/gen_rep_elems_for_cache_from_file.py
echo "Finished Task: $SLURM_ARRAY_TASK_ID"
