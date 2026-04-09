#!/bin/bash

# for future use of making an "interactive version"
export FLAGS="$1"

get_parameters_and_load_temp(){
source ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
echo "Parameters loaded; running cache generation."
mkdir $MFC_TEMP
}

setup_file_system(){
if [ ! -d "$BKPT_PATH_BASE/$NUM_BKPT" ]; then
    mkdir -p $BKPT_PATH_BASE/$NUM_BKPT
else
    echo "path: $BKPT_PATH_BASE/$NUM_BKPT exists; proceed with caution."
fi
if [ ! -d "$REP_ELEM_PATH_BASE/$NUM_BKPT" ]; then
    mkdir -p $REP_ELEM_PATH_BASE/$NUM_BKPT
else
    echo "path: $REP_ELEM_PATH_BASE/$NUM_BKPT exists; proceed with caution."
fi
}

setup_cache_run(){
echo "Running inital setup."
chmod +x ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh
#if [ "$FLAGS" == *c* ]; then 
srun --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$INITAL_TIME:00 --wait ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh $NUM_BKPT
#else
# ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh $NUM_BKPT
#fi
# Now load job info
source $MFC_TEMP/temp_job_info_for_$NUM_BKPT.sh

# check if we should keep goin'
if [ $RUN_COMPUTATION == '0' ]; then
  echo "Run aborted."
  echo "Cleaning temp files."
  rm -rf $MFC_TEMP
  exit 0
else
  echo "Run Approved - Computing $NUM_BKPT."
fi
}

setup_apptainer(){
echo "Checking for apptainer."
module load apptainer
if [ ! -f "$MFC_TARGET/cgf.sif" ]; then
  echo "Building apptainer"
  apptainer build "$MFC_TARGET/cgf.sif" "$APPTAINER_DEF_PATH"
fi
}

main(){
# get parameters and load temp
source ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/cache_gen_run_parameters.sh
echo "Parameters loaded; running cache generation."
mkdir $MFC_TEMP
for ((NUM_BKPT = $NUM_BKPT_LOWER_BOUND; NUM_BKPT <= $NUM_BKPT_UPPER_BOUND; NUM_BKPT++))
do
echo "Starting run for $NUM_BKPT."

# setup file system
if [ ! -d "$BKPTS_PATH_BASE/$NUM_BKPT" ]; then
    mkdir -p $BKPTS_PATH_BASE/$NUM_BKPT
else
    echo "path: $BKPT_PATH_BASE/$NUM_BKPT exists; proceed with caution."
fi
if [ ! -d "$REP_ELEM_PATH_BASE/$NUM_BKPT" ]; then
    mkdir -p $REP_ELEM_PATH_BASE/$NUM_BKPT
else
    echo "path: $REP_ELEM_PATH_BASE/$NUM_BKPT exists; proceed with caution."
fi

#setup_apptainer
echo "Checking for apptainer."
module load apptainer
if [ ! -f "$MFC_TARGET/cgf.sif" ]; then
  echo "Building apptainer"
  apptainer build "$MFC_TARGET/cgf.sif" "$APPTAINER_DEF_PATH"
fi

# setup_cache_run
echo "Running inital setup."
chmod +x ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh
#if [ "$FLAGS" == *c* ]; then 
sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time="$INITAL_TIME:00" --mem=$MEM --wait ./MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh $NUM_BKPT
#else
# ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/setup_cache_run.sh $NUM_BKPT
#fi
# Now load job info
source $MFC_TEMP/temp_job_info_for_$NUM_BKPT.sh

# check if we should keep goin'
if [ '$RUN_COMPUTATION' == '0' ]; then
  echo "Run aborted."
  echo "Cleaning temp files."
  rm -rf $MFC_TEMP
  exit 0
else
  echo "Run Approved - Computing function cache for  $NUM_BKPT. Disbatching jobs."
fi

# echo "Run setup complete $NUM_BKPT. Dispatching jobs."
chmod +x ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh
sbatch --array=0-$NUM_JOBS --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --time="$ALLOC_TIME_PER_JOB:00" --mem=$MEM ~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/job_function_runner.sh $NUM_BKPT

done
echo "Cleaning temp files"
rm -rf $MFC_TEMP
}

main

