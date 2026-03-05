from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import time 
import os
from numpy import std
from random import sample
import logging

inital_gen_logger = logging.getLogger(__name__)

# read bash inputs
k = int(system_args[1]) # for now keep k small, <= 9 from initial testing. 
sample_size = int(system_args[2]) # pick a value between 2 and 100;
time_per_batch = float(system_args[3]) # minutes
max_number_of_rows = int(system_args[4]) 
max_std = float(system_args[5])
try:
    which_backend = str(system_args[6]) # a string specifying backend, optional.
except IndexError:
    which_backend = None
try: 
    overhead_time = float(system_args[7]) # minutes
except:
    overhead_time = 5
try:
    run_computation_default = str(system_args[8])
except IndexError:
    run_computation_default = None
batch_time_default = 10

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
    inital_gen_logger.info("Will as about running full computation later.")

assert(k >= 2)
assert(time_per_batch >= batch_time_default)
assert(sample_size >= 2)
assert(sample_size <= 100)
assert(max_number_of_rows >= 1)
inital_gen_logger.info("Computing breakpoints...")
logging.disable()
bkpts = BreakpointComplexClassContainer(k)

# bkpts.write_data(max_rows=max_lines_in_file)

# Estimate cpu time per computation per breakpoint
num_bkpts_seqs = len(list(bkpts.get_rep_elems()))
inital_gen_logger.info(f"Number of breakpoints: {num_bkpts_seqs}.\nDetermining run time estimate.")
sample_space = sample(list(bkpts.get_rep_elems()), k = sample_size)
inital_gen_logger.info("Sample space found; computing sample points...")
start = time.time()
find_minimal_function_reps_from_bkpts(sample_space)
end = time.time()
average_sample_time = (end - start) / sample_size
inital_gen_logger.info(f"Sample computation time averages {average_sample_time:.2f} second per breakpoint sequence.")
# find time to be used per batch. 
# Infer a reasonable time estimate based on user specifications and measured data. 
# we don't want to request more time than needed on accident. 



# fix dim analysis here
estimated_cpu_time = average_sample_time * num_bkpts_seqs
inital_gen_logger.info(f"Total estimated cpu time is {estimated_cpu_time/60:.2f} minutes")
inital_gen_logger.info("Finding output size in rows and number of batches ...")
guess_number_of_batches = int((estimated_cpu_time / (time_per_batch * 60)))+1
#print(f"We will use {number_of_batches} for computation.")
number_of_rows =  int(num_bkpts_seqs/guess_number_of_batches)
if number_of_rows > max_number_of_rows:
    batch_time_less_than_estimated_computation_time = True
    number_of_rows = max_number_of_rows
else:
    batch_time_less_than_estimated_computation_time = False
number_of_batches = int(num_bkpts_seqs/number_of_rows) + 1
inital_gen_logger.info(f"Number of rows per bkpt file:  {number_of_rows}\nNumber of rows per min_fun_rep file: {number_of_rows*k} \nNumber of Batches: {number_of_batches}")
inital_gen_logger.info("Inferring batch time estimate based on parameters")
sample_std =  std(sample_space)
if overhead_time*60 > sample_std * max_std:
    time_alloc = time_per_batch + overhead_time
else:
    time_alloc = time_per_batch + (num_rows*sample_std * max_std)/60 
print(f"Allocated time batch: {time_alloc:.2f} minutes.")
print(f"Allocating {time_alloc * number_of_batches:.2f} minutes of CPU time.")

# if not set to default run or stop, ask user to decide if the program should run
if run_computation is None:
    print(f"Do you want to continue to run the computation? Y or N?")
    valid_input = False
    while not valid_input:
        cont = input()
        if cont.lower().strip(' ') == "n":
            run_computation = False
            valid_input = True
            print("Aborting run...")
        elif cont.lower().strip(' ') == "y":
            run_computation = True
            valid_input = True
            print("Run approved. Carrying on with computation.")
        else:
            print("Please input either Y or N to continue or not")
    # prompt user to set run_computation, if not already set
    

#os.putenv(RUN_JOBS, run_computation) # some type of bool  do 
if run_computation:
    # write breakpoints
    # Set ENV variables for this run
    # record file writing location and add to path.
    #os.putenv(NUM_ROWS, number_of_rows)
    #os.putenv(NUM_JOBS, number_of_batches)
    #os.putenv(ALLOC_TIME_PER_JOBS, time_alloc)
    pass
