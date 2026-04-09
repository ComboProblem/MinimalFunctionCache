from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import time
import os
from numpy import std
from random import sample
import logging

# set up logger
initial_gen_logger = logging.getLogger(__name__)
logging_level = os.getenv("LOGGING_LEVEL")
if logging_level == "debug":
    initial_gen_logger.info(f"Adjusting the default logging level to DEBUG")    
    initial_gen_logger.setLevel(logging.DEBUG)
elif logging_level == "warning":
    initial_gen_logger.info(f"Adjusting the default logging level to WARNING")
    initial_gen_logger.setLevel(logging.WARNING)
elif logging_level == "error":
    initial_gen_logger.info(f"Adjusting the default logging level to ERROR")
    initial_gen_logger.setLevel(logging.ERROR)
else:
    initial_gen_logger.setLevel(logging.INFO)
initial_gen_logger.info(f"Logging level set to {initial_gen_logger.level}")
# read parameters from cache_gen_run_parameters.sh and setup_cache_run.sh
def read_non_path_parameters(path_to_file=None):
    """
    Read and validate non path parameters from cache_gen_run_parameters.sh and setup_cache_run.sh.
    """
    k = int(os.getenv("K")) # for now keep k small, <= 9 from initial testing.
    assert(k >= 2)
    sample_size = int(os.getenv("SAMPLE_SIZE"))# pick a value between 2 and NUM_BKPTS
    assert(sample_size >= 2)
    time_per_batch = float(os.getenv("TIME_PER_BATCH")) # minutes
    assert(time_per_batch >= 1)
    max_number_of_rows = int(os.getenv("MAX_NUM_ROW")) #
    assert(max_number_of_rows >= 1)
    max_std = float(os.getenv("MAX_STD"))
    assert(max_std > 0)
    which_backend = str(os.getenv("BACKEND")) # a string specifying backend, optional.
    if which_backend == "":
        which_backend = None
    overhead_time = float(os.getenv("OVERHEAD_TIME")) # minutes
    if overhead_time == "":
        overhead_time = 5
    run_computation = os.getenv("RUN_COMPUTAION")
    if run_computation == "":
        run_computation = None
    if run_computation is not None:
        if run_computation.lower().strip(" ") == "n" or run_computation.lower().strip(" ") == "no":
            run_computation = False
        elif run_computation.lower().strip(" ") == "y" or run_computation.lower().strip(" ") == "yes":
            run_computation = True
        else:
            run_computation = None # signal to program to wait for user input.
    else:
        run_computation = None
    run_params = {"k":k, "sample_size":sample_size, "time_per_batch":time_per_batch, "max_number_of_rows":max_number_of_rows, "max_std":max_std, "which_backend":which_backend, "overhead_time":overhead_time, "run_computation":run_computation}
    initial_gen_logger.debug(f"Read Parameters: {run_params}")
    return run_params

def validate_paths(run_params):
    """
    Read path environment variables and check existence. 
    """
    paths = {"bkpts_path_base":os.getenv("BKPTS_PATH_BASE"), "rep_elem_path": os.getenv("REP_ELEM_PATH_BASE"), "temp":os.getenv("MFC_TEMP")}
    paths["bkpt_current_path"] = os.path.join(paths["bkpts_path_base"], str(run_params['k']))
    paths["bkpt_prev_path"] = os.path.join(paths["bkpts_path_base"], str(run_params['k']-1))
    parths_are_valid = True
    bad_paths = {}
    for path in paths:
        if os.path.exists(paths[path]):
            initial_gen_logger.debug(f"Path {path} with dir {paths[path]} exists.")
        else:
            initial_gen_logger.debug(f"Path {path} with dir {paths[path]} does not exists.")
            parths_are_valid = False
            bad_paths[path] = paths[path]
    if parths_are_valid is False:
        initial_gen_logger.warning(f"Some invalid path(s) provided. Listing: {bad_paths}")
    return paths
    
def gen_bkpts(paths, run_params):
    """
    Tries to use the provided specific previous breakpoints path to generate the current jobs breakpoints.
    """
    initial_gen_logger.debug("Number of files in {}\n {}".format(paths["bkpt_prev_path"], len(os.listdir(paths["bkpt_prev_path"]))))
    if len(os.listdir(paths["bkpt_prev_path"])) > 0:
        bkpts = BreakpointComplexClassContainer(run_params['k'], backend=run_params['which_backend'], manually_load_breakpoint_cache=True, folder_or_file="folder", path_to_file_or_folder=paths["bkpt_prev_path"])
        initial_gen_logger.debug("Loading previous breakpoints from {}".format(paths["bkpt_prev_path"]))
    else:
        initial_gen_logger.debug("Path {} contains no files. Generating without previous information.".format(paths["bkpt_prev_path"]))
        bkpts = BreakpointComplexClassContainer(run_params['k'], backend=run_params['which_backend'])
#    except FileNotFoundError:
#        initial_gen_logger.debug("Path {} does not exist. Generating without previous information.".format(paths["bkpt_prev_path"]))
#        bkpts = BreakpointComplexClassContainer(run_params['k'], backend=run_params['which_backend'])
    return bkpts

def estimate_time(run_params, bkpts):
    """
    Estimates the total time based on bkpts and run parameters. 
    """
    assert(run_params['sample_size'] < bkpts.num_rep_elems())
    # bkpts.write_data(max_rows=max_lines_in_file)
    # Estimate cpu time per computation per breakpoint
    initial_gen_logger.info(f"Number of breakpoints sequences: {bkpts.num_rep_elems()}.\n Determining run time estimate.")
    sample_space = sample(list(bkpts.get_rep_elems()), run_params['sample_size'])
    initial_gen_logger.info("Sample space found; computing sample points...")
    # use real time estimations
    start = time.time()
    find_minimal_function_reps_from_bkpts(sample_space, backend=run_params['which_backend'])
    end = time.time()
    average_sample_time = (end - start) / run_params['sample_size']
    initial_gen_logger.info(f"Sample computation time averages {average_sample_time:.2f} second per breakpoint sequence.")
    # find time to be used per batch.
    # Infer a reasonable time estimate based on user specifications and measured data.
    estimated_cpu_time = average_sample_time * bkpts.num_rep_elems() # total timein seconds
    initial_gen_logger.info(f"Total estimated cpu time is {estimated_cpu_time/60:.2f} minutes")
    initial_gen_logger.info("Finding output size in rows and number of batches ...")
    guess_number_of_batches = int((estimated_cpu_time / (run_params['time_per_batch'] * 60)))+1
    number_of_rows =  int(bkpts.num_rep_elems() / guess_number_of_batches)
    if number_of_rows > run_params['max_number_of_rows']:
        # batch_time_less_than_estimated_computation_time == True
        number_of_rows = run_params['max_number_of_rows']
    # else:
        # batch_time_less_than_estimated_computation_time == False
    number_of_batches = int(bkpts.num_rep_elems()/number_of_rows) + 1
    initial_gen_logger.info(f"Number of rows per bkpt file:  {number_of_rows}\nNumber of rows per min_fun_rep file: {number_of_rows * (run_params['k'] - 1)} \n Number of Batches: {number_of_batches}")
    initial_gen_logger.info("Inferring batch time estimate based on parameters")
    sample_std =  std(sample_space)
    initial_gen_logger.debug(f"Sample std: {sample_std}")
    if run_params['overhead_time']*60 > sample_std * run_params['max_std']:
        time_alloc = run_params['time_per_batch'] + run_params['overhead_time']
    else:
        time_alloc = run_params['time_per_batch'] + (number_of_rows*sample_std * run_params['max_std'])/60
    time_alloc = int(time_alloc+1)
    initial_gen_logger.info(f"Allocated time batch: {time_alloc} minutes.")
    initial_gen_logger.info(f"Allocating {time_alloc * number_of_batches:.2f} minutes of CPU time.")
    batch_info = {'number_of_rows':number_of_rows, 'number_of_batches':number_of_batches, 'time_alloc':time_alloc}
    initial_gen_logger.debug(f"Batch info: {batch_info}")
    return batch_info
# Interactive decisions to authorize the computation.
def interactive_run(run_params):
    if run_params['run_computation'] is None:
        print(f"Do you want to continue to run the computation? Y or N?")
        valid_input = False
        while not valid_input:
            cont = input()
            if cont.lower().strip(' ') == "n":
                run_params['run_computation'] = False
                valid_input = True
                initial_gen_logger.info(f"Aborting run.")
            elif cont.lower().strip(' ') == "y":
                run_params['run_computation'] = True
                valid_input = True
                initial_gen_logger.info(f"Run approved. Carrying on with computation.")
            else:
                print(f"Please input either Y or N to continue or not.")
    return run_params
# Write a results files that bash will read. 
def create_job(batch_info, bkpts, paths, run_params):
    job_info_name = "temp_job_info_for_{}.sh".format(run_params['k'])
    cwd = os.getcwd()
    if run_params['run_computation']:
        os.chdir(paths["bkpt_current_path"])
        bkpts.write_data(max_rows=run_params['max_number_of_rows'])
        os.chdir(cwd)
        os.chdir(paths["temp"])
        with open(job_info_name, "w") as run_vars_file:
            run_vars_file.write("#!/bin/bash\n")
            run_vars_file.write(f"export run_computation=1\n")
            run_vars_file.write(f"export NUM_ROWS={batch_info['number_of_rows']}\n")
            run_vars_file.write(f"export NUM_JOBS={batch_info['number_of_batches']}\n")
            run_vars_file.write(f"export ALLOC_TIME_PER_JOB={batch_info['time_alloc']}")
    else:
        os.chdir(paths["temp"])
        with open(job_info_name, "w") as run_vars_file:
            run_vars_file.write("#!/bin/bash\n")
            run_vars_file.write(f"export run_computation=0")

def __main__():
    initial_gen_logger.info("Reading Parameters and validating paths")
    run_params = read_non_path_parameters()
    paths = validate_paths(run_params)
    initial_gen_logger.info("Computing breakpoints...")
    bkpts = gen_bkpts(paths, run_params) 
    initial_gen_logger.info("Estimating batch information")
    batch_info = estimate_time(run_params, bkpts)
    initial_gen_logger.debug("Checking interactive run.")
    run_params = interactive_run(run_params)
    initial_gen_logger.info("Writing Job Parameters...")
    create_job(batch_info, bkpts, paths, run_params)
    initial_gen_logger.info("Finished.")

__main__()
