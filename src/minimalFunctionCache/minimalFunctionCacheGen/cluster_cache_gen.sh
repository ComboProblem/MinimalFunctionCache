#!/bin/bash -x

# Get run parameters.
source ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh

echo "Parameters loaded; running cache generation."
mkdir $MFC_TEMP
for NUM_BKPT in {$NUM_BKPT_LOWER_BOUND..$NUM_BKPT_UPPER_BOUND}
do
echo "Starting run for $NUM_BKPT."
# Set up file system; should add some logic to 
mkdir $MFC_TARGET
mkdir $MFC_TARGET/Breakpoints
mkdir $MFC_TARGET/RepElems
mkdir $MFC_TARGET/Breakpoints/$NUM_BKPT
mkdir $MFC_TARGET/RepElems/$NUM_BKPT


# load module(s)
# module purge
echo "Checking for apptainer."
module load apptainer
if [ ! -f cgf.sif ]; then
  echo ""
  apptainer build cgf.sif $APPTAINER_DEF_PATH
fi

# Run the inital set up.
echo "Running inital setup."
chmod +x ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh
~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh

sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$INITAL_TIME:00  ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh

# Now load job info
source $MFC_TEMP/temp_job_info_for_$NUM_BKPT.sh

# check if we should keep goin'
if [$RUN_COMPUTATION == '0']; then
  echo "Run aborted."
  echo "Cleaning temp files."
  rm -rf $MFC_TEMP
  exit 0
else
  echo "Starting run for $NUM_BKPT"
fi

chmod +x ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh
sbatch --array=0-$NUM_JOBS --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --time=$ALLOC_TIME_PER_JOB:00 --mem=$MEM ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh

done
echo "Cleaning temp files"
rm -rf $MFC_TEMP


