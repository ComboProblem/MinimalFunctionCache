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

export BKPTS_PATH = "~/MinimalFunctionCache/Breakpoints/$NUM_BKPT"
mkdir ~/MinimalFunctionCache/TEMP
export MFC_TEMP = "~/MinimalFunctionCache/TEMP"
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
e
  # specify what path(s) to bind inside the apptainer
  export $APPTAINER_BIND=$PWD

  # check for an apptainer(singularity) image (.sif) and build/pull one if it doesn't exist
  # Assumes that Apptainer.def is in PWD.

  apptainer run cgf.sif cache_job_setup.py NUM_BKPT SAMPLE_SIZE TIME_PER_BATCH MAX_STD NUM_ROW BACKEND

# verify
echo NUM_ROWS: $NUM_ROWS
echo NUM_JOBS: $NUM_JOBS
echo ALLOC_TIME_PER_JOB: $ALLOC_TIME_PER_JOB
  # check new env vars are set

  # add a continue option/flag always yes, always no, ask when finished (default behavior).

  # Here is where we make a loop for excute jobs. We
