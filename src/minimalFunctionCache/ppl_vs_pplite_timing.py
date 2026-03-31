from cutgeneratingfunctionology.igp import *
from minimalFunctionCache import utils

def backend_timing(num_bkpt, backend=None, number_of_sf=6, max_value_poly_clock_time=None):
    """
    Times value polyhedron calcuations.
    """
    if max_value_poly_clock_time is None:
        max_value_poly_clock_time = 2**256-1 # time in second, humans and the universe will be dead by then 
    bkpts = minimal_function_cache_loader(num_bkpt, "breakpoints")
    timing_data = []
    for bkpt in bkpts:
        for f_index in range(1, num_bkpt):
            start = time.time()
            value_nnc_polyhedron_value_cords(list(bkpt), f_index, backend)
            end = time.time()
            data = round(start-end, number_of_sf)
            timing_data.append(data)
    return data


def run_timing_tests(clock_time_per_test=1000):
    cache_info = minimal_function_cache_info()
    ppl_data = {}
    pplite_data = {}
    for num_bkpt in cache_info["avail_bkpts"]:
        ppl_data[num_bkpt] = backend_timing(num_data, max_value_poly_clock_time=clock_time_per_test)
        pplite_data[num_bkpt] = backend_timing(num_data, backend="pplite", max_value_poly_clock_time=clock_time_per_test)
    
    
