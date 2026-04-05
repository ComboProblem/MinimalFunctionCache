#!/bin/bash -x

#SLURM JOB INFORMATION
export CLUSTER_ACCOUNT=math-grp
export PARTITION=high
export MEM=4gb
export INITAL_TIME=120 # in minutes
# Set cache computation parameters
export NUM_BKPT_LOWER_BOUND=6 # integers, lower <= upper bound
export NUM_BKPT_UPPER_BOUND=10 
export SAVE_BKPTS=true  # bash bool
export SAVE_REP_ELEMS=true # bash bool
export SUBMIT_TO_GITHUB=true # bash bool
# For estimating cpu time.
export SAMPLE_SIZE=10 # increasing sample size increases the initial run time which could long with current implementation.
export TIME_PER_BATCH=60 # measured in minutes
export OVERHEAD_TIME_PER_BATCH=5 # in minutes
export MAX_STD=5 # used in time estimages
export MAX_NUM_ROW=1000 # in each breakpoint file.
export BACKEND="pplite" #pplite or None
export OVERHEAD_TIME=5 # in minutes
export LOGGING_LEVEL="" #debug, warning, or error, anything else will default to info
# Default Run Option
export RUN_COMPUTAION=yes # after breakpoints computation dispatch to the cluster
# Locations of paths
export MFC_TEMP="~/MinimalFunctionCache/TEMP"
export MFC_TARGET="~/MinimalFunctionCache/src/minimalFunctionCache"
export BKPTS_PATH="$MFC_TARGET/Breakpoints/$NUM_BKPT"
export REP_ELEM_PATH="$MFC_TARGET/RepElems/$NUM_BKPT"
export APPTAINER_DEF_PATH="~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/Apptainer.def"

