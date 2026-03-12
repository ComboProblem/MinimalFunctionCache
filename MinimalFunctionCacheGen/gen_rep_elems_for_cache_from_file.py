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
rep_elm_gen_logger.info(f"Starting computations for {job_number}.")
try:
    PiMin_worker = PiMinContContainer(k, load_bkpt_data=input_file_name, backend=backend)
    rep_elm_gen_logger.info(f"Computations from {job_number} has finished. Writing data.")
    output_file_name = f"Pi_Min_{k}_part_{job_number}"
    os.chdir(os.getenv("REP_ELEM_PATH"))
    PiMin_worker.write(output_file_name)
except FileNotFoundError:
    rep_elm_gen_logger.warning(f"{input_file_name} not found in {os.getenv("BKPTS_PATH")}")
