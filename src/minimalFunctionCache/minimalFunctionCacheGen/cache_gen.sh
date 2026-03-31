#!/bin/bash -x
#-x is useful for debugging.

# Get run parameters.
source ~/cache_gen_run_parameters.sh

# Set up temp folders.
mkdir ~/MinimalFunctionCache/TEMP
mkdir ~/MinimalFunctionCache/TEMP/Breakpoints
mkdir ~/MinimalFunctionCache/TEMP/Breakpoints/$NUM_BKPT

if [ $SAVE_BKPTS ]; then
  if [ ! -d ~/MinimalFunctionCache/MinimalFunctionCache/Breakpoints/$NUM_BKPTS ]; then
    mkdir ~/MinimalFunctionCache/MinimalFunctionCache/Breakpoints/$NUM_BKPTS
  fi
fi


# load module(s)
# module purge
module load apptainer
if [ ! -f cgf.sif ]; then
  apptainer build cgf.sif ~/Apptainer.def
fi

# Run the inital set up.
# chmod +x ~/setup_cache_run.sh
# sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$  apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/cache_job_setup.py $NUM_BKPT $SAMPLE_SIZE $TIME_PER_BATCH $MAX_NUM_ROW $MAX_STD $BACKEND $OVERHEAD_TIME_PER_BATCH $RUN_COMPUATION
sbatch --partition=low --ntasks=1 --cpus-per-task=1 --time=$  apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/cache_job_setup.py $NUM_BKPT $SAMPLE_SIZE $TIME_PER_BATCH $MAX_NUM_ROW $MAX_STD $BACKEND $OVERHEAD_TIME_PER_BATCH $RUN_COMPUATION

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


mkdir ~/MinimalFunctionCache/TEMP/RepElems
mkdir ~/MinimalFunctionCache/TEMP/RepElems/$NUM_BKPT
# if [ $SAVE_BKPTS ]; then
#   if [ ! -d ~/MinimalFunctionCache/MinimalFunctionCache/Breakpoints/$NUM_BKPTS ]; then
#     mkdir ~/MinimalFunctionCache/MinimalFunctionCache/Breakpoints/$NUM_BKPTS
#   fi
# fi
chmod +x ~/job_function_runner.sh

apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/gen_rep_elems_for_cache_from_file.py
