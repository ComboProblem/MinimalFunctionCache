#!/bin/bash

source ./MinimalFunctionCache/src/minimalFunctionCacheGen/cache_gen_run_parameters.sh

module load apptainer
apptainer exec "$MFC_TEMP/cgf_mem.sif" python -m memory_profiler ./MinimalFunctionCache/src/minimalFunctionCacheGen/MemoryProfiling/profile_breakpoints.py
