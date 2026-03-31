#!/bin/bash -x
source ~/cache_gen_run_parameters.sh
source ~/MinimalFunctionCache/TEMP/temp_job_info.sh

#SBATCH --array=0-$NUM_JOBS
#SBATCH --cpus-per-task=1
#SBATCH --time=$ALLOC_TIME_PER_JOB:00
#SBATCH --mem=4GB
#SBATCH --partition=low

module load apptainer
apptainer run cgf.sif python3 ~/MinimalFunctionCache/MinimalFunctionCacheGen/gen_rep_elems_for_cache_from_file.py
