from cutgeneratingfunctionology.igp import *
from minimalFunctionCache.utils import *
from memory_profiler import LogFile
import sys
# create logger
logger = logging.getLogger('memory_profile_log')
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.FileHandler("MinimalFunctionCache/src/minimalFunctionCacheGen/MemoryProfiling/four_bkpt_profile.log")
fh.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(fh)
sys.stdout = LogFile('memory_profile_log', reportIncrementFlag=False)

six_bkpts = minimal_function_cache_loader(6, "breakpoints")

# This version removed logging code and documentation; for examining memory profile of code actually important. 
@profile    
def make_rep_bkpts_with_len_n_profile_edition(n, k=1, bkpts=None, backend=None):
    new_bkpts = []
    if n < 2:
        raise ValueError("n>=2")
    if k == n and bkpts is not None:
        return bkpts
    if k == n and bkpts is None:
        raise ValueError("k<n")
    if k == 1 and bkpts is None:
        bkpts = [[0]]
    for bkpt in bkpts:
        new_bkpts += add_breakpoints_and_find_equiv_classes(nnc_poly_from_bkpt_sequence(bkpt, backend=backend).upstairs())
    new_bkpts = unique_list(new_bkpts)
    k += 1
    if k == n:
        return new_bkpts
    else:
        return make_rep_bkpts_with_len_n(n, k, new_bkpts, backend)
        
seven_bkpts = make_rep_bkpts_with_len_n_profile_edition(7, 6, six_bkpts, "pplite")
