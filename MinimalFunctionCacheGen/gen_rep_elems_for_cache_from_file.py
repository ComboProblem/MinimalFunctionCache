from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import os
import logging
rep_elm_gen_logger = logging.getLogger(__name__)
rep_elm_gen_logger.setLevel(logging.INFO)

k = int(system_args[1])
job_number = os.
path = os.path.curdir
rep_elm_gen_logger.info(f"Hello from job {job_number}.")
input_file_name = None
output_file_name = None
# here we should check that the file exists.
if input_file_name is in path:
    rep_elm_gen_logger.info(f"Loading {input_file_name}...")
    rep_elm_gen_logger.info(f"Starting computations for {job_number}.")
    PiMin_worker = PiMinContContainer(load_bkpt_data=input_file_name)
    rep_elm_gen_logger.info(f"Computations from {job_number} has finished. Writing data.")
    PiMin_worker.write(output_file_name)
else:
    rep_elm_gen_logger.error(f"{input_file_name} is not found in {path}.}")
    
