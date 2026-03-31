#!/bin/bash -x

#SLURM JOB INFORMATION
# export CLUSTER_ACCOUNT=<account_name>
# export PARTITION=<partition>

# Set cache computation parameters
export NUM_BKPT=6
export SAVE_BKPTS=true
# For estimating cput time.
export SAMPLE_SIZE=10 # increasing sample size increases the inital run time which could lenghtly wiht current implementation.
export TIME_PER_BATCH=60 # mesured in minutes
export OVERHEAD_TIME_PER_BATCH=5 # in minutes
export MAX_STD=5 # used in time estimages
export MAX_NUM_ROW=1000 # in each breakpoint file.
export BACKEND="pplite"

# Locations of paths
export MFC_TEMP="~/MinimalFunctionCache/TEMP"
export BKPTS_PATH="$MFC_TEMP/Breakpoints/$NUM_BKPT"
export REP_ELEM_PATH="$MFC_TEMP/RepElems/$NUM_BKPT"
