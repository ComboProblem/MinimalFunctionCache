#!/bin/bash

setup_apptainer(){
echo "Checking for apptainer."
module load apptainer
if [ ! -f "MinimalFunctionCache/src/minimalFunctionCacheGen/cgf.sif" ]; then
  echo "Building apptainer"
  apptainer build "MinimalFunctionCache/src/minimalFunctionCacheGen/cgf.sif" "MinimalFunctionCache/src/minimalFunctionCacheGen/Apptainer.def"
else
   echo "Cutgeneratingfunctiology image already exists."
fi
}

main(){

source ./MinimalFunctionCache/src/minimalFunctionCacheGen/cache_gen_run_parameters.sh
echo "Parameters loaded; running cache generation."
mkdir $MFC_TEMP
#setup_apptainer
echo "Checking for apptainer."
module load apptainer
if [ ! -f "MinimalFunctionCache/src/minimalFunctionCacheGen/cgf.sif" ]; then
  echo "Building apptainer"
  apptainer build "$MFC_TEMP/cgf_mem.sif" "MinimalFunctionCache/src/minimalFunctionCacheGen/MemoryProfilingApptainer_prof.def"
else
   echo "cutgeneratingfunctiology image already exists."
fi

chmod +x ~/MinimalFunctionCache/src/minimalFunctionCacheGen/MemoryProfiling/memory_profile_run.sh
sjob --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --time="720:00" --mem=50G ~/MinimalFunctionCache/src/minimalFunctionCacheGen/memory_profile_run.sh

rm -rf $MFC_TEMP
}

main

