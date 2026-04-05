#!/bin/bash -x

echo "Loading Cluster ENV" 
source ~/MinimalFunctionCache/src/minimalFunctionCache/parametericCutGen/cluster_enviroment.sh

mkdir EXP_TEMP

echo "Checking for apptainer."

module load apptainer
if [ ! -f optimal_cut.sif ]; then
  echo ""
  apptainer build optimal_cut.sif $APPTAINER_DEF_PAT

echo "Configuring experiments"
sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=SETUP_TIME:00 exp_setup.sh

echo "Dispatching experiments."

