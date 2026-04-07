import os
import csv
import logging
import importlib.resources as importlib_resources
from sage.rings.rational_field import QQ


utils_logger = logging.getLogger(__name__)

def minimal_function_cache_dir():
    """
    Returns directory information for accessing the minimal function cache.
    """
    function_cache_path = importlib_resources.files("minimalFunctionCache")
    breakpoint_base_path = function_cache_path / 'Breakpoints'
    rep_elem_base_path = function_cache_path / 'RepElems'
    return {"function_cache_path":function_cache_path, "breakpoint_base_path":breakpoint_base_path, "rep_elem_base_path":rep_elem_base_path }

def minimal_function_cache_info():
    """
    Retuns the known minimal funciton cache information. 
    """
    fun_cache_dir = minimal_function_cache_dir()
    minimal_function_cache_info = {"avail_bkpts": [] , "avail_rep_elems": [] }
    for entry in fun_cache_dir["breakpoint_base_path"].iterdir():
        minimal_function_cache_info["avail_bkpts"].append(int(entry.name))
    for entry in fun_cache_dir["rep_elem_base_path"].iterdir():
        minimal_function_cache_info["avail_rep_elems"].append(int(entry.name))
    return minimal_function_cache_info

def minimal_function_cache_loader(n, breakpoints_or_rep_elems, prototype=QQ):
    """
    Loads minimal function cache elements as a list with each repereseative element values having the type of prototype, by default, QQ.

    Prototype is a callable that that takes on input a string with format "a/b" where a and b are integers.
    Prototype can return anything.
    
    EXAMPLES::
    
    >>> from MinimalFunctionCache.cacheUtils import cache_loader
    >>> bkpts_for_2 = cache_loader(2, "breakpoints") # not tested
    >>> len(bkpts_for_2) # not tested
    3
    >>> rep_elems_for_2 = cache_loader(2, "rep_elems") # not tested
    >>> len(rep_elems_for_2) # not tested
    3
    >>> def python_float_prototype(rational_as_string):
                a, b = rational_as_string.split("/")
                return float(a)/float(b)
    >>> rep_elems_for_2_float = cache_loader(2, "rep_elems", python_float_prototype) # not tested
    >>> all([isinstace(elem, float) for elem in rep_elems_for_2_float[0])) # not tested
    True
    """
    fun_cache_dir = minimal_function_cache_dir()
    cache_info = minimal_function_cache_info()
    if breakpoints_or_rep_elems.strip(" ").lower() == "breakpoints":
        if n not in cache_info["avail_bkpts"]:
            utils_logger.info(f"Breakpoint cache {n} not found")
            raise ValueError(f"The breakpoint cache for {n} has not been computed.")
        bkpts = []
        bkpt_path = fun_cache_dir["breakpoint_base_path"] / str(n)
        for file in list(bkpt_path.glob('*.csv')):
            with open(file, newline='') as csvfile:
                file_reader = csv.reader(csvfile)
                for row in file_reader:
                    bkpts.append([prototype(data) for data in row])
        return bkpts
    elif breakpoints_or_rep_elems.strip(" ").lower() == "rep_elems":
        if n not in cache_info["avail_rep_elems"]:
            utils_logger.info(f"Rep elem cache {n} not found")
            raise ValueError(f"The function cache for {n} has not been computed.")
        rep_elems = []
        rep_elem_path = fun_cache_dir["rep_elem_base_path"] / str(n)
        for file in list(rep_elem_path.glob('*.csv')):
            with open(file, newline='') as csvfile:
                file_reader = csv.reader(csvfile)
                for row in file_reader:
                    bkpt = [QQ(data) for data in row[0].strip("[]").split(",")]
                    val = [QQ(data) for data in row[1].strip("[]").split(",")]
                    rep_elems.append((bkpt, val))
        return rep_elems
    else:
        raise ValueError("A cache has not been loaded, check spelling and inputs.")




