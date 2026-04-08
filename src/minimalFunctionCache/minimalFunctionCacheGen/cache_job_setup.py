from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import time
import os
from numpy import std
from random import sample
import logging

# set up logger
inital_gen_logger = logging.getLogger(__name__)
logging_level = os.getenv("LOGGING_LEVEL")
if logging_level == "debug":
    inital_gen_logger.INFO(f"Adjusting the default logging level to DEBUG")
    inital_gen_logger.setLevel(logging.DEBUG)
elif logging_level == "warning":
    inital_gen_logger.INFO(f"Adjusting the default logging level to WARNING")
    inital_gen_logger.setLevel(logging.WARNING)
elif logging_level == "error":
    inital_gen_logger.INFO(f"Adjusting the default logging level to ERROR")
    inital_gen_logger.setLevel(logging.ERROR)
else:
    inital_gen_logger.setLevel(logging.INFO)
# read bash inputs
# TODO: Rewrite as reading enviroment variables.
k = int(system_args[1]) # for now keep k small, <= 9 from initial testing.
sample_size = int(system_args[2]) # pick a value between 2 and NUM_BKPTS
time_per_batch = float(system_args[3]) # minutes
max_number_of_rows = int(system_args[4]) #
max_std = float(system_args[5])
try:
    which_backend = str(system_args[6]) # a string specifying backend, optional.
except IndexError:
    which_backend = None
try:
    overhead_time = float(system_args[7]) # minutes
except IndexError:
    overhead_time = 5
try:
    run_computation_default = str(system_args[8])
except IndexError:
    run_computation_default = None
if inital_gen_logger.level is logging.DEBUG:
    logging.debug({"k":k, "sample_size":sample_size, "time_per_batch":time_per_batch, "max_number_of_rows":max_number_of_rows, "max_std":max_std, "which_backend":which_backend, "overhead_time":overhead_time, "run_computation_default":run_computation_default})
inital_gen_logger.info("Setting up SLURM job parameters.")
if run_computation_default is not None:
    if run_computation_default.lower().strip(" ") == "n" or run_computation_default.lower().strip(" ") == "no":
        run_computation = False
        inital_gen_logger.info("The full computation will not run.")
    elif run_computation_default.lower().strip(" ") == "y" or run_computation_default.lower().strip(" ") == "yes":
        run_computation = True
        inital_gen_logger.info("The full computation will run after inferring job parameters.")
    else:
        inital_gen_logger.info("No valid argument for default choice to run full computation or not. Will ask again later.")
        run_computation = None # signal to program to wait for user input.
else:
    run_computation = None
    inital_gen_logger.info("Will ask about running full computation later.")
assert(k >= 2)
assert(time_per_batch >= 1)
assert(sample_size >= 2)
assert(max_number_of_rows >= 1)
inital_gen_logger.info("Computing breakpoints...")
os.environ["BKP_PRV_DIR"]="./MinimalFunctionCache/src/minimalFuncitonCache/Breakpoints/{k-1}"
try:
    if len(os.listdir(os.getenv('BKP_PRV_DIR'))) > 0:
        bkpts = BreakpointComplexClassContainer(k, backend=which_backend, manually_load_breakpoint_cache=True, file_or_foler="folder", path_to_file_or_foler=os.getenv('BKP_PRV_DIR'))
        inital_gen_logger.debug(f"Loading previuos breakpoints from {os.getenv('BKP_PRV_DIR')}")
    else:
        inital_gen_logger.debug(f"Path {os.getenv('BKP_PRV_DIR')} contains no files. Generating breakpoints without prior infromation.")
        bkpts = BreakpointComplexClassContainer(k, backend=which_backend)
except FileNotFoundError:
    inital_gen_logger.debug(f"Path {os.getenv('BKP_PRV_DIR')} does not exist. Generating breakpoints without prior infromation.")
    bkpts = BreakpointComplexClassContainer(k, backend=which_backend)

assert(sample_size < bkpts.num_rep_elems())
# bkpts.write_data(max_rows=max_lines_in_file)
# Estimate cpu time per computation per breakpoint
inital_gen_logger.info(f"Number of breakpoints sequences: {bkpts.num_rep_elems()}.\nDetermining run time estimate.")
sample_space = sample(list(bkpts.get_rep_elems()), k = sample_size)
inital_gen_logger.info("Sample space found; computing sample points...")
# use real time estimations
start = time.time()
find_minimal_function_reps_from_bkpts(sample_space, backend=which_backend)
end = time.time()
average_sample_time = (end - start) / sample_size
inital_gen_logger.info(f"Sample computation time averages {average_sample_time:.2f} second per breakpoint sequence.")
# find time to be used per batch.
# Infer a reasonable time estimate based on user specifications and measured data.
estimated_cpu_time = average_sample_time * num_bkpts_seqs # total timein seconds
inital_gen_logger.info(f"Total estimated cpu time is {estimated_cpu_time/60:.2f} minutes")
inital_gen_logger.info("Finding output size in rows and number of batches ...")
guess_number_of_batches = int((estimated_cpu_time / (time_per_batch * 60)))+1
number_of_rows =  int(num_bkpts_seqs/guess_number_of_batches)
if number_of_rows > max_number_of_rows:
    # batch_time_less_than_estimated_computation_time == True
    number_of_rows = max_number_of_rows
# else:
    # batch_time_less_than_estimated_computation_time == False
number_of_batches = int(num_bkpts_seqs/number_of_rows) + 1
inital_gen_logger.info(f"Number of rows per bkpt file:  {number_of_rows}\nNumber of rows per min_fun_rep file: {number_of_rows*k} \nNumber of Batches: {number_of_batches}")
inital_gen_logger.info("Inferring batch time estimate based on parameters")
sample_std =  std(sample_space)
if overhead_time*60 > sample_std * max_std:
    time_alloc = time_per_batch + overhead_time
else:
    time_alloc = time_per_batch + (number_of_rows*sample_std * max_std)/60
time_alloc = int(time_alloc+1)
inital_gen_logger.info(f"Allocated time batch: {time_alloc} minutes.")
inital_gen_logger.info(f"Allocating {time_alloc * number_of_batches:.2f} minutes of CPU time.")
# Interactive decisions to authorize the computation.
# TODO: Make this read a bash var for a -i flag when running the script.
if run_computation is None:
    print(f"Do you want to continue to run the computation? Y or N?")
    valid_input = False
    while not valid_input:
        cont = input()
        if cont.lower().strip(' ') == "n":
            run_computation = False
            valid_input = True
            print(f"Aborting run...")
        elif cont.lower().strip(' ') == "y":
            run_computation = True
            valid_input = True
            print(f"Run approved. Carrying on with computation.")
        else:
            print(f"Please input either Y or N to continue or not.")
# Write a results files that bash will read. 
job_info_name = "temp_job_info_for_{}.sh".format(k)
if run_computation:
    os.chdir(os.getenv("BKPTS_PATH"))
    bkpts.write_data(max_rows=max_number_of_rows)
    os.chdir(os.getenv("MFC_TEMP"))
    with open(job_info_name, "w") as run_vars_file:
        run_vars_file.write("#!/bin/bash\n")
        run_vars_file.write(f"export RUN_COMPUTATION=1\n")
        run_vars_file.write(f"export NUM_ROWS={number_of_rows}\n")
        run_vars_file.write(f"export NUM_JOBS={number_of_batches}\n")
        run_vars_file.write(f"export ALLOC_TIME_PER_JOB={time_alloc}")
else:
    os.chdir(os.getenv("MFC_TEMP"))
    with open(job_info_name, "w") as run_vars_file:
        run_vars_file.write("#!/bin/bash\n")
        run_vars_file.write(f"export RUN_COMPUTATION=0")
os.exit()
