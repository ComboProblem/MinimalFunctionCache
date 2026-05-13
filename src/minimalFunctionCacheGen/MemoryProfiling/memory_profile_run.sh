#!/bin/bash

source ./MinimalFunctionCache/src/minimalFunctionCacheGen/cache_gen_run_parameters.sh

module load apptainer
apptainer run "$MFC_TEMP/cgf_mem.sif" python3 memory_profiler profile_breakpoints.py
