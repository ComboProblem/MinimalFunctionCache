import os
import csv
import importlib.resources as importlib_resources
from sage.rings.rational_field import QQ


# Set up paths
_function_cach_path = importlib_resources.files("cacheUtils")
_breakpoint_base_path = _function_cach_path / 'Breakpoints'
_rep_elem_base_path = _function_cach_path / 'RepElems'

cache_info = {"avail_bkpts": [] , "avail_rep_elems": [] }

# lets parse the cache.
for entry in _breakpoint_base_path.iterdir():
    cache_info["avail_bkpts"].append(int(entry.name))

for entry in _rep_elem_base_path.iterdir():
    cache_info["avail_rep_elems"].append(int(entry.name))

def cache_loader(n, breakpoints_or_rep_elems, prototype=QQ):
    """
    Loads cache elements as a list with each repereseative element values having the type of prototype, by default, QQ

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
    if breakpoints_or_rep_elems.strip(" ").lower() == "breakpoints":
        if n not in cache_info["avail_rep_elems"]:
            raise ValueError(f"The breakpoint cache for {n} has not been computed.")
        bkpts = []
        bkpt_path = _breakpoint_base_path / str(n)
        for file in list(bkpt_path.glob('*.csv')):
            with open(file, newline='') as csvfile:
                file_reader = csv.reader(csvfile)
                for row in file_reader:
                    bkpts.append([prototype(data) for data in row])
        return bkpts
    elif breakpoints_or_rep_elems.strip(" ").lower() == "rep_elems":
        if n not in cache_info["avail_bkpts"]:
            raise ValueError(f"The function cache for {n} has not been computed.")
        rep_elems = []
        rep_elem_path =  _rep_elem_base_path / str(n)
        for file in list(rep_elem_path.glob('*.csv')):
            with open(file, newline='') as csvfile:
                file_reader = csv.reader(csvfile)
                file_reader = csv.reader(csvfile)
                for row in file_reader:
                    bkpts.append([prototype(data) for data in row])
        return rep_elems
    else:
        raise ValueError("A cache has not been loaded, check spelling and inputs")




