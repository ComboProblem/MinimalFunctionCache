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
    def __init__(self, MIP_objective):
        pass
    def __call__(self, cut):
        pass


class cutGenerationDomain:
    r"""PiMin<=k, possibly with paramaterized constraints.
    constraints is an interable of BSA objects
    get_cells is the indeded call method. Should allow users to 
    """
    def __init__(self, k, possible_bkpts_as_params=[], constraints=[]):
        self._num_bkpts = k
        self._possible_bkpts_as_params = possible_bkpts_as_params # should this be a copy?
        # If no possible breakpoints are given as parameters, read this as no constraints on the domain. 
        # if some possible breakpoints are given, we assume these as possible constraints on the space, but not as necessary constraints on the space
        # necessary constraints should be given via add constraint. 
        if k > max_bkpts and len(self._possible_bkpts_as_params) == 0: 
            raise ValueError("The minimal function cache requested has not been computed. Please compute the minimal funciton cache for {}.".format(k))
        elif k <= max_bkpts:
            self._PiMin = PiMinContContainer(self._num_bkpts, load_rep_elem_data=sys_info.path_to_rep_elem_data(k)) # sys_info here has a method to collect the file location and file names and passes this into PiMinConstinaterTo Load
        else: # k > max_bkpts and bkpts_as_params
            self._PiMin = None
            # Generate slices of cells based on possible_bkpts_as_params
        for con in constraints:
            if not isinstance(con, BasicSemialgebraicSet_base):
                raise ValueError("Constraints should should be represented as BSAs")
                # should also check that maps are compatiable or can be make assinged into the correct field. There is some detail checking to be done.
        self._constraints = constraints
    def __repr__(self):
        return "The Space of Con. Minimal Functions with at most {} breakpoints and constrained by {}".format
    # def is_cell_description(self):
        # pass
    # def is_manifold_description(self):
        # pass
        
    def get_cells(self):
        r"""iterate over all cells (either explicilty defined or infered) with intersected constraints. Outputs a BSA.
        """
        if self._PiMin is not None:
            yield self._PiMin.get_semialgebraic_sets().intersection(self._constraints)
        else:
            for bkpt in Subsets(self._possible_bkpts_as_params, self._num_bkpts):
                for f_index in range(1, self._num_bkpts)
                    yield generate_cell_slice_from_bkpt_params(bkpt, f_index).intersection(self._constraints)
    # def get_manifold(self): #To do, once I write manifold methods. 
        # raise NotImplementedError
    # def topology(self, cell):
        # pass #open, closed, semi open? Is this needed? 
    # def paramaterize_point(self, bkpt, val):
        # return piecewise_function_from_breakpoints_and_values(bkpt+[1], val+[0])
    def is_element_of_domain(self, pwl_or_bkpt_and_val):
        pass
    def add_constraint(self, constraint):
        self._constraints.append(constraint) # probs should keep as a BSA
    def remove_constrain(self, constraint):
        pass
    def add_possible_breakpoint_value(self, b):
        self._possible_bkpts_as_params.append(b)
    def generate_cell_slice_from_bkpt_params(self, bkpt, f_index):
        assert(len(bkpt) ==  self._num_bkpts)
        return value_nnc_polyhedron(bkpt, f_index) #write a verison which includes makes a consistent number of parameters. 


class cutGenerationSolverBase:
    r"""
    A base class for interfacing with solvers. cutGenerationDomains can be expressed a finite collection of semialgebraic cells. 
    This has methods translate semialgebraic cell or manifold descriptions of data from the cutGenerationDomain to some external solver.
    The solve and __init__ methods are written to be compatiable with the cut generation problem  and provide an optimal solution to the problem 
    provided that the solver and data conversion methods have been written. 
    """
    def __init__(self, cutGenerationDomain, cutScore):
        pass
    def solve(self, **options):
        pass
    def write_linear_constraints_from_bsa_for_solver(self, bsa): # think about aspects of exactness; 
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError
    def write_nonlinear_constraints_from_bsa(self, bsa):
        r"""
        Given a BSA with non linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError
    # def write_manifold_constraints(self, manifold):
        # pass
    def solver_linear_solve(self, constraints, objective, min_or_max,  **options):
        r"""
        Interface to solver's solve method. Should output an optimal value and point to the problem min/max c^Tx s.t. Ax <= b. 
        """
        raise NotImplementedError
    def solver_nonlinear_solve(self, constraints, objective, **options):
        raise NotImplementedError
    def is_solver_exact(self):
        r"""
        Return True if solver is does exact rational computations, false otherwise. 
        """
        pass

class cutGenerationProblem:
    r"""
    Define the problem min cutScore(pi, MIP) s.t. pi in PiMin <=k
    """
    def __init__(self, max_bkpts, cut_scoring_method, MIP, cut_domain, solver, **solving_parameters):
        self._status = None
        if isinstance(cut_domain, cutGenerationDomain):
            self._domain = cut_domain # TODO: add plain text alias for domains types
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
        




