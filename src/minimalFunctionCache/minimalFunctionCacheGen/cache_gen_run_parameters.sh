#!/bin/bash

#SLURM JOB INFORMATION
export CLUSTER_ACCOUNT=math-grp
export PARTITION=high
export MEM=1gb
declare -i INITAL_TIME=5 # in minutes
# Set cache computation parameters
declare -i NUM_BKPT_LOWER_BOUND=2 # integers, lower <= upper bound
declare -i NUM_BKPT_UPPER_BOUND=2
export SAVE_REP_ELEMS=true # bash bool
export SUBMIT_TO_GITHUB=true # bash bool
# For estimating cpu time.
export SAMPLE_SIZE=2 # increasing sample size increases the initial run time which could long with current implementation.
export TIME_PER_BATCH=10 # measured in minutes
export OVERHEAD_TIME_PER_BATCH=1 # in minutes
export MAX_STD=3 # used in time estimates
export MAX_NUM_ROW=1000 # in each breakpoint file.
export BACKEND="pplite" #pplite or None
export OVERHEAD_TIME=5 # in minutes
export LOGGING_LEVEL="debug" #debug, warning, or error, anything else will default to info
# Default Run Option
export RUN_COMPUTAION=yes # after breakpoints computation dispatch to the cluster
# Locations of paths
export MFC_TEMP="./MinimalFunctionCache/TEMP"
export MFC_TARGET="./MinimalFunctionCache/src/minimalFunctionCache"
export BKPTS_PATH_BASE="$MFC_TARGET/Breakpoints"
export REP_ELEM_PATH_BASE="$MFC_TARGET/RepElems"
export APPTAINER_DEF_PATH="./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/Apptainer.def"
