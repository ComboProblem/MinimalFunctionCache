from cutgeneratingfunctionology.igp import *
from minimalfunctioncache import sys_info
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
# sys_info is an object that contains knows how to read data from the minimal function cache.

# current defaults
max_bkpts = sys_info.max_bkpts
solver_mode_default = "full"
row_processing_order_default = "lex"
cut_scoring_method_default = "parallelism"

class UnsetData(Exception):
    pass

class ParamaterizationError(Exception):
    pass

class abstractCutScore:
    r"""
    Class factory for cut optimization objective functions, aka cut scoring metrics. 
    """
    @classmethod
    def __init__(cls, **kwrds):
        pass

    @abstractmethod
    def cut_score(cls, cut, field=None):
        r"""
        Assume cut is a list like object of elements of the given field (by default QQ). 
        Cut is assumed to be the intersection cut generated from some minimal function. 
        In particular, we are assuming that has the form \sum_{j\in N} pi(bar(a_ij))^Tx_j 
        = \sum_{j\in N} pi_p(bar(a_ij))x_j \geq 1 = \pi_p(bar(b_i))= pi_p(lambda_findex) = 1.
        I.e. (0, pi)(x_B,x_N) \geq 1.
        
        Output of cut_score should be an element of the specified field (or should allow sage to coherse the elements). 
        """
        raise NotImplementedError

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

    @classmethod
    def cut_obj_type(cls)
         cls._cut_obj_type

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
            raise TypeError("name a valid cut score or ")

    def __init__(self, cut_score, **kwrds):
        r"""
        Initialize the cut scoring function. 
        
        MIP data, type conversion method, and objective type, can be specified in initialization or modified
        with with get and set methods. 
        """
        self._cut_score = cut_score(**kwrds)
        super().__init__()
        self._MIP_objective = None
        self._MIP_row = None
        self._sage_to_solver_type = None
        self._obj_type = None
    
    def __call__(self, parameters):
        r"""
        parameters is a list like object with even length of at most 2k.
        parameters = (bkpt, val) represents a parameterize element of Pimin<=k, i.e. a point in RR^{2n} such that
        the first n elements are the breakpoints, and the remainder of the values of the list are the values
        of the parameterized function, pi_(bkpt,val) \in Pimin<=k.
        
        It is assumed parameters is satisfied the assumptions stated. No errors will be raised. 
        """
        # parameters is and element of B\times V \se R^2n
        # we assume parameters are of the solver type.
        # we preform the composition of maps which takes us from the parameter space to
        # minimal funciton space. 
        # Internally, we use exact rational arithmetic as that is what works for the exact
        # constructions of parametric functions.
        # then convert the output to the correct solver type finishing the composition of 
        # maps. 
        if self._MIP_objective is None:
            raise UnsetData("Set MIP_objective before use of CutScore.")
        if self._MIP_row is None:
            raise UnsetData("Set MIP_row before use of CutScore.")
        if self._sage_to_solver_type is None:
            raise UnsetData("Set sage_to_solver_type before use of CutScore.")
        bkpt = [QQ(b) for b in parameters[:len(parameters)/2]]
        val = [QQ(v) for v in parameters[len(parameters)/2:]]
        pi = piecewiese_function_from_bkpts_and_values(bkpt + [1], val + [0])
        row_data =  self._cut_score.get_MIP_row()
        cut = [pi(QQ(bar_a_ij)) for bar_a_ij in row_data]
        result = self.get_sage_to_solver_type()(self._cut_score.cut_score(cut))
        return result    

class Parallelism(abstractCutScore):
    """
    Normalized cut parallelism score.
    """
    def cut_score(cls, cut, field=None):
        obj_norm = vector(cls._MIP_objective).norm()
        cut_norm = vector(cut).norm()
        dot_product = vector(cls._MIP_objective).row()*vector(cut).column()
        return dot_product[0]/(obj_norm*cut_norm)
        
class SteepestDirection(abstractCutScore)
    """
    Steepest direction score. 
    """
    def cut_score(cls, cut, field=None):
        dot_product = vector(cls._MIP_objective).row()*vector(cut).column()
        return dot_product[0]
        

class cutGenerationDomain:
    r"""PiMin<=k, possibly with paramaterized constraints.
   
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
        for bkpt in Subsets(self.possible_bkpt_params):
            for f_index in range(1, self._num_bkpts)
                yield generate_cell_slice_from_bkpt_params(bkpt, f_index)
    
    def specify_f(self, f):
        """
        Assumes all elements in the cutGenerationDomain satisfy the condition that pi(f) = 1.
        """
        self._f = f 
    
    def current_f(self):
        return self._f

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
        self._constraints.pop(constraint)

    def add_possible_breakpoint_value(self, b):
        self._possible_bkpt_params.append(b)

    def generate_cell_slice_from_bkpt_params(self, bkpt, f_index):
        assert(len(bkpt) ==  self._num_bkpts)
        return value_nnc_polyhedron(bkpt, f_index) #write a version which includes makes a consistent number of parameters. 


class cutGenerationSolverBase:
    r"""
    A base class for interfacing with solvers. cutGenerationDomains can be expressed a finite collection of semialgebraic cells. 
    This has methods translate semialgebraic cell or manifold descriptions of data from the cutGenerationDomain to some external solver.
    The solve and __init__ methods are written to be compatiable with the cut generation problem  and provide an optimal solution to the problem 
    provided that the solver and data conversion methods have been written. 
    """
    def __init__(self, cut_gen_domain, cut_score, cut_problem_options, **solver_options):
        if not isinstance(cut_gen_domain, cutGenerationDomain):
            raise ValueError("cut_gen_domain is required to be an instance of cutGenerationDomain")
        if not isinstance(cut_score, cutScoreBase):
            raise ValueError("cut_score is requried to be an instance of cutScoreBase")
        self._cut_gen_domain = cut_gen_domain
        self._cut_score = cut_score
        self._cut_score.set_sage_to_solver_type(self.sage_to_solver_type)

    def solve(self, MIP, min_or_max, **solver_options):
        r"""Solves the paramaterized problem 
        
        min/max cutScore s.t. 
        options are options to be passed into the solver. 
        """
        # for now assuming problem is in the form of min 
        current_domain_obj_val = inf
        domain_problem_solution = None
        
        # ADD Options to which row get solved in the MIP.relaxed_basis()
        for fractional_row in MIP.relaxed_basis():
            self._cut_score.set_MIP_row(fractional_row)
            f_constraint = "fractional_row[b_i] == f"
            val_constraint =  "pi(f) == 1"
            self._cut_gen_domain.specify_f(fractional_row[b_i]) #replace with actually BSA constraint.

### solve problem for a set row. 
        for subdomain in self._cut_gen_domain.get_cells():
            objective_fun = self._cut_score # the call method here is a function from R^{2n} (parameter space) to R
            if subdomain.is_linear() #reprlace with correct BSA command:
                subdomain_solver_constraints = self.write_linear_constraints_from_bsa_for_solver(subdomain)
                subdomain_problem_val, subdomain_problem_solution = self.solver_linear_solve(subdomain_solver_constraints, objective_fun, min_or_max, **solver_options)
            else: # non-linear case
                subdomain_solver_constraints = self.write_nonlinear_constraints_from_bsa_for_solver(subdomain)
                subdomain_problem_val, subdomain_problem_solution = self.solver_nonlinear_solve(subdomain_solver_constraints, objective_fun, min_or_max, **solver_options)                
            if subdomain_problem_val < current_domain_obj_val:
                current_domain_obj_val = subdomain_problem_val
                domain_problem_solution = subdomain_problem_solution
        result = {}
        result["objective_value"] = current_domain_obj_val
        result["parameterized_solution"] = domain_problem_solution
        return result

    @staticmethod
    def write_linear_constraints_from_bsa_for_solver(bsa): # think about aspects of exactness; 
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError

    @staticmethod
    def write_nonlinear_constraints_from_bsa(self, bsa):
        r"""
        Given a BSA with non linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError

    # def write_manifold_constraints(self, manifold):
        # pass

    @staticmethod
    def solver_linear_solve(self, constraints, objective, min_or_max,  **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b. 
        
        Should return optimal objective value, optimal objective solution.
        """
        raise NotImplementedError

    @staticmehtod
    def solver_nonlinear_solve(self, constraints, objective, min_or_max, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polynomial and at least 1 p_i has degree larger than 1.  
        """
        raise NotImplementedError


    @staticmethod
    def sage_to_solver_type(self, sage_field_element, field=None):
        raise NotImplementedError
    

class cutGenerationProblem:
    r"""
    MIP has data c^Tx st. Ax<=b; x>=0; A is n times m; x is len m.  
    Define the problem min cutScore(pi, MIP) s.t. pi in PiMin <=k. 
    
    Option: row selection: all_rows, lex rows, random row, subset (find optimal cut over fractional row(s))
    Option: algorithm = full space; bkpt as param: full (all combinatorial data, if needed k<n-m),  max, lex (selects lexicographically first parameter data k< n-m), or rand (randomly select parameter data)
    Option: num_bkpt = full space: k <= max_bkpts; bkpt as param: full, lex rand, k <= n-m; max, k = n-m
    Option: cut_score = parallelism, steepestdirection, scip, or custom
    Option: cut_gen_solver = scipy or custom
    Option: multithread = not_implemented
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
