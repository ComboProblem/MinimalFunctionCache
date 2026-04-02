#!/bin/bash -x

#SLURM JOB INFORMATION
# export CLUSTER_ACCOUNT=<account_name>
export PARTITION=high
export MEM=4gb
export INITAL_TIME=120
# Set cache computation parameters
export NUM_BKPT=4
export SAVE_BKPTS=true
export SAVE_REP_ELEMS=true
export SUBMIT_TO_GITHUB=true
# For estimating cpu time.
export SAMPLE_SIZE=10 # increasing sample size increases the inital run time which could lenghtly wiht current implementation.
export TIME_PER_BATCH=60 # mesured in minutes
export OVERHEAD_TIME_PER_BATCH=5 # in minutes
export MAX_STD=5 # used in time estimages
export MAX_NUM_ROW=1000 # in each breakpoint file.
export BACKEND="pplite"
export LOGGING_LEVEL="" #debug, warning, or error, anything else will default to info

# Locations of paths
export MFC_TEMP="~/MinimalFunctionCache/TEMP"
export BKPTS_PATH="$MFC_TEMP/Breakpoints/$NUM_BKPT"
export REP_ELEM_PATH="$MFC_TEMP/RepElems/$NUM_BKPT"
