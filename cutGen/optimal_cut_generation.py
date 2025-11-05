from cutgeneratingfunctionology.igp import *
from minimalfunctioncache import sys_info
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
# sys_info is an object that contains knows how to read data from the minimal function cache.


# current defaults
max_bkpts = sys_info.max_bkpts
solver_mode_default = "full"
row_processing_order_default = "lex"
cut_scoring_method_default = "scip_default"


class abstractCutScore:
    r"""
    Class factory for cut optimization objective functions, aka cut scoring metrics. 
    """
    @classmethod
    def __init__(cls, name, MIP_objective, MIP_row, sage_to_solver_type):
        cls._name = name
        cls._MIP_objective = MIP_objective
        cls._MIP_row = MIP_row
        cls._sage_to_solver_type =  sage_to_solver_type

    @abstractmethod
    def cut_score(cls, cut):
        r"""
        Assume cut is a list of sage rationals. 
        
        The call function here is to evaluate the current parameterized cut to in the cutGenerationSolver.  
        The cut as the form \sum_{j\in N} pi(bar(a_ij))^Tx_j= \sum_{j\in N} pi_p(bar(a_ij))x_j \geq 1 = \pi_p(bar(b_i))= pi_p(lambda_findex) = 1.
        We have the  constraint lambda_findex = b_i on pi Min assumed to be holding at this point. 
        B corrosponds to the basis in a current LP relaxation basis of MIP.
        N is the non basic variables of the LP relaxation. 
        This is where the cut scoring method (which is assumed to be fixed relative to a fixed objective function of the MIP)
        The particular cut scoring method should overwrite the call function and replace it with assuming that cut is a list like object 
        representing the cut. cutGenSolver is responsible for translating types to solver compaitable types. 
        """
        raise NotImplementedError
    
    @classmethod
    def name(cls):
        return cls._name

    @classmethod
    def set_MIP_objective(cls, new_objective):
        cls._MIP_objective = new_objective

    @classmethod
    def get_MIP_objective(cls):
        return cls._MIP_objective

    @classmethod
    def get_MIP_row(cls):
        return cls._MIP_row

    @classmethod
    def set_MIP_row(cls, new_row):
        cls._MIP_row = new_row

    @classmethod
    def get_sage_to_solver_type(cls):
        return cls._sage_to_solver_type

    @classmethod
    def set_sage_to_solver_type(cls, new_conversion):
        cls.__sage_to_solver_type = new_conversion

class CutScore:
    @staticmethod
    def __classcall__(cls, name=None, **kwrds):
        r"""

        """
        if name == "parallelism" or name is None:
            return super().__classcall__(cls, cut_score=Parallelism)
        if issubclass(name, abstractCutScore):
            return super().__classcall__(cls, cut_score=name, **kwrds)
        else:
            raise TypeError("BOO")

    def __init__(self, cut_score, **kwrds):
        self._cut_score = cut_score(**kwrds)
        super().__init__()
    
    def __call__(self, parameters):
        # parameters is and element of B\times V \se R^2n
        # we assume parameters are of the solver type.
        # we preform the composition of maps which takes us from the parameter space to
        # minimal funciton space. 
        # Internally, we use exact rational arithmetic as that is what works for the exact
        # constructions of parametric functions.
        # then convert the output to the correct solver type finishing the composition of 
        # maps. 
        bkpt = [QQ(b) for b in parameters[:len(parameters)/2]]
        val = [QQ(v) for v in parameters[len(parameters)/2:]]
        pi = piecewiese_function_from_bkpts_and_values(bkpt + [1], val + [0])
        row_data =  self._cut_score.get_MIP_row()
        cut = [pi(QQ(bar_a_ij)) for bar_a_ij in row_data]
        result = self.get_sage_to_solver_type()(self._cut_score.cut_score(cut))
        return result    


class cutGenerationDomain:
    r"""PiMin<=k, possibly with paramaterized constraints.
   
    constraints is an interable of BSA objects
    get_cells is the indeded call method. Should allow users to 
    """
    def __init__(self, k, possible_bkpt_params=[], constraints=[]):
        self._num_bkpts = k
        self._possible_bkpt_params = breakpointSequence(possible_bkpt_params)
        # If no possible breakpoints are given as parameters, read this as no constraints on the domain. 
        # if some possible breakpoints are given, we assume these as possible constraints on the space, but not as necessary constraints on the space
        # necessary constraints should be given via add constraint. 
        if k > max_bkpts and len(self.possible_bkpt_params) == 0: 
            raise ValueError("The minimal function cache requested has not been computed. Please compute the minimal funciton cache for {}.".format(k))
        elif k <= max_bkpts:
            self._PiMin = PiMinContContainer(self._num_bkpts, load_rep_elem_data=sys_info.path_to_rep_elem_data(k)) # sys_info here has a method to collect the file location and file names and passes this into PiMinConstinaterTo Load
        else: # k > max_bkpts and bkpts_as_params
            self._PiMin = None
        for con in constraints:
            if not isinstance(con, BasicSemialgebraicSet_base):
                raise ValueError("Constraints should should be represented as BSAs")
                # should also check that maps are compatiable or can be make assinged into the correct field. There is some detail checking to be done.
        # Think about how BSA constraints sould be written and type normalizaton. 
        self._constraints = constraints

    def __repr__(self):
        return "The Space of Con. Minimal Functions with at most {} breakpoints and constrained by {}.".format(self._num_bkpts, self._constraints)
    # def is_cell_description(self):
        # pass
    # def is_manifold_description(self):
        # pass
        
    def get_cells(self):
        r"""iterate over all cells (either explicilty defined or infered) with intersected constraints. Outputs a BSA.
        """
        if self._PiMin is not None:
            cells = self._PiMin.get_semialgebraic_sets()
        else:
            cell = None
        for bkpt in Subsets(self.possible_bkpt_params, self._num_bkpts):
            for f_index in range(1, self._num_bkpts)
                yield generate_cell_slice_from_bkpt_params(bkpt, f_index)

    # def get_manifold(self): #To do, once I write manifold methods. 
        # raise NotImplementedError

    # def topology(self, cell):
        # pass #open, closed, semi open? Is this needed? 

    # def paramaterize_point(self, bkpt, val):
        # return piecewise_function_from_breakpoints_and_values(bkpt+[1], val+[0])

    def is_element_of_domain(self, pwl_or_bkpt_and_val):
        pass

    def add_constraint(self, constraint):
        self._constraints.append(constraint) 

    def remove_constrain(self, constraint):
        pass

    def add_possible_breakpoint_value(self, b):
        self._possible_bkpt_params.append(b)

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
    def __init__(self, cut_gen_domain, cut_score)):
        if not isinstance(cut_gen_domain, cutGenerationDomain):
            raise ValueError("cut_gen_domain is required to be an instance of cutGenerationDomain")
        if not isinstance(cut_score, cutScoreBase):
            raise ValueError("cut_score is requried to be an instance of cutScoreBase")
        self._cut_gen_domain = cut_gen_domain
        self._cut_score = cut_score
        

    def solve(self, MIP, **options):
        r"""Solves the paramaterized problem options are options to be passed into the solver. 
        """
        for subdomain in self._cut_gen_domain.get_cells():
            if subdomain.is_linear() #reprlace with correct BSA command:
                subdomain_solver_constraints = self.write_linear_constraints_from_bsa_for_solver(subdomain)
                objective_fun = self._cut_score # the call method here is a function from R^{2n} (parameter space) to R
        pass

    @staticmethod
    def write_linear_constraints_from_bsa_for_solver(bsa): # think about aspects of exactness; 
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
        Interface to solver's min/max f(x) s.t. Ax<=b. 
        """
        raise NotImplementedError

    def solver_nonlinear_solve(self, constraints, objective, **options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polymomial and at least 1 p_i has degree larger than 1.  
        """
        raise NotImplementedError
    def sage_to_solver_type(self, cut):
        raise NotImplementedError
    
    

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
        

# 100 percent using SCIP
