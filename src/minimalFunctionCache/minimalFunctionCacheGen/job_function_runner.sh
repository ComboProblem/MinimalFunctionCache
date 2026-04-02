#!/bin/bash -x

source ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
source ~/MinimalFunctionCache/TEMP/temp_job_info.sh

module load apptainer
apptainer run cgf.sif python3 ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/gen_rep_elems_for_cache_from_file.py
