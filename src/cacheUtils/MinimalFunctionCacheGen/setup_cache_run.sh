#!/bin/bash -x
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=4GB
#SBATCH --partition=low

echo "Running on $(hostname)"
echo "Setting up run."
source ~/cache_gen_run_parameters.sh
module load apptainer
apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/cache_job_setup.py $NUM_BKPT $SAMPLE_SIZE $TIME_PER_BATCH $MAX_NUM_ROW $MAX_STD $BACKEND
