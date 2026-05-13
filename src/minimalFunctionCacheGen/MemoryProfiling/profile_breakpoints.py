from cutgeneratingfuctionology.igp import *
from minimalFunctionCache import utils

six_bkpts = minimal_function_cache_loader(6, "breakpoints")

out_file =  open("./src/minimalFunctionCacheGen/MemoryProfiling/seven_bkpt_profile.log", "w+")

@profile(stream=out_file)
def make_rep_bkpts_with_len_n(n, k=1, bkpts=None, backend=None):
    r"""
    Produce representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n.

    INPUT:
    - n, integer, maximum length of breakpoint sequence.
    - k, assumed length of breakpoint sequences in ``bkpts``.
    - bkpts, list of breakpoint sequences of length k.

    OUTPUT: A list of representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n extrapolated from bkpts.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: make_rep_bkpts_with_len_n(2)
    [(0, 1/2), (0, 13/18), (0, 5/18)]

    The number of representative elements grows quickly::

    sage: bkpts_rep_with_len_3 = make_rep_bkpts_with_len_n(3)
    sage: len(bkpts_rep_with_len_3)
    34

    Previous computations can be reused::

    sage: bkpts_rep_with_len_4 = make_rep_bkpts_with_len_n(4, 3, bkpts_rep_with_len_3)
    sage: len(bkpts_rep_with_len_4)
    329
    """
    # Matthias has suggested looking at a directed tree.
    # An alternative approach would be to look into using a (graded) lattice as a data structure.
    # We have bkpt \leq bkpt' if and only if dim(NNC(bkpt)) \leq dim(NNC(bkpt'))
    # and embed(NNC(bkpt), dim(NNC(bkpt')) \leq NNC(bkpt') in the poset of NNC polyhedra.
    # This might speed up/less the load of verifying unquiness of cells which is the time bounding task here.
    new_bkpts = []
    if n < 2:
        raise ValueError("n>=2")
    if k == n and bkpts is not None:
        minimal_funciton_cell_description_logger.warning(f"Initial inputs suggest that the bkpts provided are already correct. Returning the initial bkpts.")
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
        minimal_funciton_cell_description_logger.info(f"Breakpoints of length {n} have been generated.")
        return new_bkpts
    else:
        minimal_funciton_cell_description_logger.info(f"Breakpoints of length {k} have been generated. Now generating breakpoints of length {k+1}.")
        return make_rep_bkpts_with_len_n(n, k, new_bkpts, backend)
        
seven_bkpts = make_rep_bkpts_with_len_n(7, 6, six_bkpts, "pplite")
