#!/bin/bash -x

#SLURM JOB INFORMATION
export CLUSTER_ACCOUNT=math-grp
export PARTITION=high
export MEM=4gb
export SETUP_TIME=10 # in minutes
export MAX_EXP_TIME=60  # in minutes
export OVERHEAD_TIME=5 # in minutes
#PATHS
export APPTAINER_DEF_PATH="~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/Apptainer.def"
export DATA_OUT_FILE_PATH=""
export EXP_TEMP="~/MinimalFunctionCache/TEMP"
#Experimental Parameters
export EXP_PARAM_SCIP_PATH="~/MinimalFunctionCache/src/parametricCutGen/paramFiles/scip_experimental_settings.toml"
export EXP_PARAM_SCIPY_PATH="~/MinimalFunctionCache/src/parametricCutGen/paramFiles/SciPy_experimental_settings.toml"
export EXP_PARAM_PATH="~/MinimalFunctionCache/src/parametricCutGen/paramFiles/experimental_parameters.toml"
