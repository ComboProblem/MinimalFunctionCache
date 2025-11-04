from cutgeneratingfunctionology.igp import *
from minimalfunctioncache import sys_info
# sys_info is an object that contains knows how to read data from the minimal function cache.


# current defaults
max_bkpts = sys_info.max_bkpts
solver_mode_default = "full"
row_processing_order_default = "lex"
cut_scoring_method_default = "scip_default"

class cutScoreBase:
    r"""
    Abstract base class for cut optimization objective functions, aka cut scoring metrics. 
    """
    pass


class cutGenerationDomain:
    r""" A subset of Pimin<=k, possibly with constraints.
    """
    def __init__(self, k, **constraints):
        pass

    def is_domain_linear(self):
        pass


class cutGenerationSolverBase:
    r"""
    Abstract base class for interfacing with linear and non linear solvers from different packages. 
    """
    def __init__(self):
        pass
    def get_linear_solver(self):
        pass
    def get_nonlinear_solver(self):
        pass

class cutGenerationProblem:
    r"""
    Define the problem min cutScore(pi, MIP) s.t. pi in cutGeneratingFunctionDomain. 
    """
    def __init__(self, max_bkpts, cut_scoring_method, MIP, cut_domain, solver, **solving_parameters):
        self._status = None
        if isinstance(cut_domain, cutGenerationDomain):
            self._domain = cut_domain # TODO: add plain text alias for domains
        if self._domain.is_linear():
            self._solver = solver.get_linear_solver()
        else:
            self._solver = solver.get_nonlinear_solver()
        self._cut_scoring_method = cutScore(cut_scoring_method)
        # TODO: Process solver_options here, for now we write defaults.
        self._solve_mode = solver_mode_default # "full", "single row", 
        self._row_processing_order = "lex"
        self._verbocity = 0 # 0, none, 1 - minimal, 2 - full. 
        self._MIP = MIP
    
    def __repr__(self):
        pass
    
    def solve(self):
        pass
        




