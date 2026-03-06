#!/bin/bash

# This script is assumed to be ran on the host computer and can send slurm jobs.

# Set cache computation parametesr
NUM_BKPT = 2
# Estimate CPU time required.
SAMPLE_SIZE = 2 # increasing sample size increases the inital run time which could lenghtly wiht current implementation.
TIME_PER_BATCH = 60 # mesured in minutes
OVERHEAD_TIME_PER_BATCH = 5 # in minutes
MAX_STD = 5 # used in time estimages
MAX_NUM_ROW = 1000 # in each breakpoint file.
BACKEND = "pplite"

mkdir ~/MinimalFunctionCache/TEMP
export MFC_TEMP = "~/MinimalFunctionCache/TEMP"
export BKPTS_PATH = "$MFC_TEMP/Breakpoints/$NUM_BKPT"

  # SBATCH --
  # SBATCH --partition=low
  # SBATCH --ntasks=1
  # SBATCH --cpus-per-task=1
  # SBATCH --mem=8G
  # SBATCH --gpus=0
  # SBATCH --job-name=cgf_apptainer_run_test
  # SBATCH --time=5:00

  # load module(s)
  module purge
  module load apptainer

  # specify what path(s) to bind inside the apptainer
  export $APPTAINER_BIND=$PWD
# typically home directory is shared over the nextwork
# file systems are realitive slow
  # check for an apptainer(singularity) image (.sif) and build/pull one if it doesn't exist
  # Assumes that Apptainer.def is in PWD.

  apptainer run cgf.sif cache_job_setup.py NUM_BKPT SAMPLE_SIZE TIME_PER_BATCH MAX_STD NUM_ROW BACKEND

# load set up variables
$MFC_TEMP/temp_job_info
# check if we should keep goin'
if [$RUN_COMPUTATION==0]; then
  echo "Run aborted."
  exit 0
fi
