from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import time 
import os
from random import sample


# read bash inputs
k = int(system_args[1]) # for now keep k small, <= 9 from initial testing. 
sample_size = int(system_args[2]) # pick a value between 2 and 100; 
time_per_batch = int(system_args[3]) # minutes
max_number_of_rows = int(system_args[4])

max_std = int(system_args[5])
which_backend = system_args[6]# a string specifying backend, optional. 
cpu_overhead_time_per_batch = 5 # overhead time allowed for each batch. 
batch_time_default = 10



assert(k >= 2)
assert(time_per_batch >= batch_time_default)
assert(sample_size >= 2)
assert(sample_size <= 100)
assert(max_number_of_rows >= 1)

logging.disable()
bkpts = BreakpointComplexClassContainer(k)

# bkpts.write_data(max_rows=max_lines_in_file)

# Estimate cpu time per computation per breakpoint
num_bkpts_seqs = len(list(bkpts.get_rep_elems()))
print(f"Number of breakpoints: {num_bkpts_seqs}.\\ Determining run time estimate.")
sample_space = sample(list(bkpts.get_rep_elems()), k = sample_size)
print("Sample space found; computing sample points")
start = time.time()
find_minimal_function_reps_from_bkpts(sample_space)
end = time.time()
average_sample_time = (end - start) / sample_size
print(f"Sample computation time averages {average_sample_time:.2f} second per breakpoint sequence.")
# find time to be used per batch. 
# Infer a reasonable time estimate based on user specifications and measured data. 
# we don't want to request more time than needed on accident. 



# fix dim analysis here
estimated_cpu_time = average_sample_time * num_bkpts_seqs
print(f"Total estimated cpu time is {estimated_cpu_time/60:2f} minutes")
print("Finding output size in rows and number of batches ...")
guess_number_of_batches = int((estimated_cpu_time / (time_per_batch * 60)))+1
#print(f"We will use {number_of_batches} for computation.")
 =  int(num_bkpts_seqs/guess_number_of_batches)
if number_of_rows > max_number_of_rows:
    batch_time_less_than_estimated_computation_time = True
    number_of_rows = max_number_of_rows
else:
    batch_time_less_than_estimated_computation_time = False
number_of_batchs = int(num_bkpts_seqs/number_of_rows) + 1
print(Number of rows per bkpt file:  {number_of_rows}\\ Number of rows per min_fun_rep file: {number_of_rows*n} \\ Number of Batches: {number_of_batches}")
print("Inferring batch time estimate based on parameters")
# compute sample standard deviation
if overhead_time*60 > sample_std * max_std:
    time_alloc = time_per_batch + overhead_time 
else:
    time_alloc = time_per_batch + (num_rows*sample_std * max_std)/60 # each batch is assumed to run in at most time_per_batch. We cover error up to max_std*sample_std in run time per ba
print(f"Allocated time batch: {time_alloc} minutes")
print(f"Allocating {time_alloc * number_of_batches} minutes of CPU time.")



# Set ENV variables for this run.
os.putenv(NUM_ROWS, number_of_rows)
os.putenv(NUM_BATCH, number_of_batches)
os.putenv(ALLOC_TIME_PER_BATCH, time_alloc)

# lets get the system location information here

# okay lets write to the system where the jobs can get files from the host. 
