# from sys import argv as system_args
from cutgeneratingfunctionology.igp import *
import os
import logging
rep_elm_worker_logger = logging.getLogger(__name__)
logging_level = os.getenv("LOGGING_LEVEL")
if logging_level == "debug":
    rep_elm_worker_logger.info(f"Adjusting the default logging level to DEBUG")    
    rep_elm_worker_logger.setLevel(logging.DEBUG)
elif logging_level == "warning":
    rep_elm_worker_logger.info(f"Adjusting the default logging level to WARNING")
    rep_elm_worker_logger.setLevel(logging.WARNING)
elif logging_level == "error":
    rep_elm_worker_logger.info(f"Adjusting the default logging level to ERROR")
    rep_elm_worker_logger.setLevel(logging.ERROR)
else:
    rep_elm_worker_logger.setLevel(logging.INFO)
rep_elm_worker_logger.info(f"Logging level set to {rep_elm_worker_logger.level}")

def set_up_job_params():
    job_number = int(os.getenv("SLURM_ARRAY_TASK_ID"))
    backend =  str(os.getenv("BACKEND"))
    k = int(os.getenv("K"))
    input_file_name = f"bkpts_of_len_{k}_part_{job_number}.csv"
    job_params = {"job_number": job_number, "backend": backend, "k":k, "input_file_name":input_file_name}
    rep_elm_worker_logger.debug(f"{job_params}")
    return job_params

def setup_and_validate_paths(job_params):
    paths = {"bkpts_path_base":os.getenv("BKPTS_PATH_BASE"), "rep_elem_path_base": os.getenv("REP_ELEM_PATH_BASE"), "temp":os.getenv("MFC_TEMP")}
    paths["bkpt_path"] = os.path.join(paths["bkpts_path_base"], str(job_params['k']))
    paths["rep_elem_path"] = os.path.join(paths["rep_elem_path_base"], str(job_params['k']))
    paths["cwd"] = os.getcwd()
    paths_are_valid = True
    bad_paths = {}
    for path in paths:
        if os.path.exists(paths[path]):
            rep_elm_worker_logger.debug(f"Path {path} with dir {paths[path]} exists.")
        else:
            rep_elm_worker_logger.debug(f"Path {path} with dir {paths[path]} does not exists.")
            paths_are_valid = False
            bad_paths[path] = paths[path]
    if paths_are_valid is False:
        rep_elm_worker_logger.error(f"All paths should exists. We have the following bad paths: {bad_paths}")
    return paths

def dispatch_worker(job_params, paths):
    try:
        PiMin_worker = PiMinContContainer(job_params["k"], manually_load_function_cache=True, file_or_folder="file", path_to_file_or_folder=os.path.join(paths["bkpt_path"], job_params["input_file_name"]), breakpoints_or_rep_elems="breakpoints")
        rep_elm_worker_logger.info(f"Computations from {job_params["job_number"]} has finished. Writing data.")
    except FileNotFoundError:
        rep_elm_worker_logger.warning(f"{job_params["input_file_name"]} not found in {paths["bkpt_path"]}")
    output_file_name = f"Pi_Min_{job_params["k"]}_part_{job_params["job_number"]}"
    os.chdir(paths["rep_elem_path"])
    PiMin_worker.write_data(output_file_name)

def __main__():
    job_params = set_up_job_params()
    paths = setup_and_validate_paths(job_params)
    rep_elm_worker_logger.info(f"Dispatching worker for job {job_params["job_number"]}.")
    dispatch_worker(job_params, paths)

__main__()
