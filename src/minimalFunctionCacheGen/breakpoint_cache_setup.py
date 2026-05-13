from cutgeneratingfunctionorgy.igp import *
import os
import time
import tomllib

"""
Estimates for distrbuted computations from previous breakpoints computations. 
"""
def validate_breakpoint_paths():
    """
    Read and validate existance of paths.
    
    Returns: dictonary of paths.
    """
    paths = {}
    return paths

def read_breakpoint_params(paths):
    """
    Read current cache computaiton parameters.
    
    Input: paths - dictionary of current path information
    Ouput: params - dictonary of read parameters.
    """
    return params

def write_breakpoint_metadata(paths, params):
    """
    Write metadata for breakpoint cache as a toml file. 

    Metadata is max length of each file, number of files, number of breakpoints.
    """
    pass

def read_breakpoint_metadata(path_to_metadata):
    """
    Output: dictionary of metadata information; using 
    """
    breakpoint_metadata = {}
    return breakpoint_metadata

def write_inital_block_generation_job_script(paths, params, previous_breakpoint_metadata):
    """
    Write a bash script to distach cluster jobs to distrubte generation of the next breakpoints.
    From the previous breakpoint files are treaded as "blocks" and are assigned to a row/column in square like
    "block" grid. The job script works on generation of and joinining of row of the underlying grid.
    Rows are made of indepenet blocks so this generation based on the row blocks can be distrubted.
    
    This writes to the temp folder. Results are targeted TEMP/NUM_BKPT
    """
    block_row_length = int(round(math.sqrt(previous_breakpoint_metadata["number_of_files"])))


def write_row_block_collection_job_script(paths, params, previous_breakpoint_metadata, block_row_length):
    """
    Write a bash script to dispatch cluster jobs which works on joining block rows. 

    To ensure uniqueness every row is checked against every other. 

    Jobs are assigned based on triangular numbers. We have (number_of_block_rows-1)number_of_block_rows/2 to check to ensure uniqueness.
    At most 2*number_of_bloc_rows batchs with number_of_block_rows/2 jobs are required. 
    With this we want to ensure a block_row can only be read/write in one job per batch (to prevent I/O waiting).


    Script is written to temp folder.
    """
    pass

def create_write_row_files_job_script(paths, params):
    """
    
    """


