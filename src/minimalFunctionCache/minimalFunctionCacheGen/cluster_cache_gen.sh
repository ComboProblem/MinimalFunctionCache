#!/bin/bash -x
#-x is useful for debugging.

# Get run parameters.
source ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh

# Set up file system
mkdir ~/MinimalFunctionCache/TEMP
mkdir ~/MinimalFunctionCache/TEMP/Breakpoints
mkdir ~/MinimalFunctionCache/TEMP/Breakpoints/$NUM_BKPT
mkdir ~/MinimalFunctionCache/TEMP/RepElems
mkdir ~/MinimalFunctionCache/TEMP/RepElems/$NUM_BKPT


# load module(s)
# module purge
module load apptainer
if [ ! -f cgf.sif ]; then
  apptainer build cgf.sif ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/Apptainer.def
fi

# Run the inital set up.
chmod +x ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh
~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh
# sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$  apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/cache_job_setup.py $NUM_BKPT $SAMPLE_SIZE $TIME_PER_BATCH $MAX_NUM_ROW $MAX_STD $BACKEND $OVERHEAD_TIME_PER_BATCH $RUN_COMPUATION
#sbatch --partition=$PARTITION --ntasks=1 --cpus-per-task=1 --time=$INITAL_TIME:00  ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh

ls ~/MinimalFunctionCache/TEMP
# Now load job info
source ~/MinimalFunctionCache/TEMP/temp_job_info.sh

# check if we should keep goin'
if [$RUN_COMPUTATION == '0']; then
  echo "Run aborted."
  echo "Cleaning temp files."
  rm -rf ~/MinimalFunctionCache/TEMP
  exit 0
else
  echo "Starting Run."
fi


chmod +x ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh
sbatch --array=0-$NUM_JOBS --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --time=$ALLOC_TIME_PER_JOB:00 --mem=$MEM ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh

if [ $SAVE_BKPTS ]; then
  if [ ! -d ~/MinimalFunctionCache/src/minimalFunctionCache/Breakpoints/$NUM_BKPTS ]; then
    mkdir ~/MinimalFunctionCache/src/minimalFunctionCache/Breakpoints/$NUM_BKPTS
  fi
  rm ~/MinimalFunctionCache/src/minimalFunctionCache/Breakpoints/$NUM_BKPTS/*.csv
  mv ~/MinimalFunctionCache/TEMP/Breakpoints/$NUM_BKPT  ~/MinimalFunctionCache/src/minimalFunctionCache/Breakpoints/$NUM_BKPTS
fi

if [ $SAVE_REP_ELEMS ]; then
  if [ ! -d ~/MinimalFunctionCache/src/minimalFunctionCache/RepElems/$NUM_BKPTS ]; then
    mkdir ~/MinimalFunctionCache/src/minimalFunctionCache/RepElems/$NUM_BKPTS
  fi
  rm ~/MinimalFunctionCache/src/minimalFunctionCache/RepElems/$NUM_BKPTS/*.csv
  mv ~/MinimalFunctionCache/TEMP/RepElems/$NUM_BKPTS ~/MinimalFunctionCache/src/minimalFunctionCache/RepElems/$NUM_BKPTS
fi


rm -rf ~/MinimalFunctionCache/TEMP
    
