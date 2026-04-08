# from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import os
import logging
rep_elm_gen_logger = logging.getLogger(__name__)
rep_elm_gen_logger.setLevel(logging.INFO)
job_number = int(os.getenv("SLURM_ARRAY_TASK_ID"))
backend =  str(os.getenv("BACKEND"))
k = int(os.getenv("NUM_BKPT"))
input_file_name = f"bkpts_of_len_{k}_part_{job_number}.csv"
os.chdir(os.getenv("BKPTS_PATH"))
try:
    rep_elm_gen_logger.info(f"Starting computations for job {job_number}.")
    PiMin_worker = PiMinContContainer(k, manually_load_function_cache=True, file_or_folder="file", path_to_file_or_folder=input_file_name, breakpoints_or_rep_elems="breakpoints")
    rep_elm_gen_logger.info(f"Computations from {job_number} has finished. Writing data.")
except FileNotFoundError:
    rep_elm_gen_logger.warning(f"{input_file_name} not found in {os.getenv("BKPTS_PATH")}")
output_file_name = f"Pi_Min_{k}_part_{job_number}"
os.chdir(os.getenv("REP_ELEM_PATH"))
PiMin_worker.write_data(output_file_name)
