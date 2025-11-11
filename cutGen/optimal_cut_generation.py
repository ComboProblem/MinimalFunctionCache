from cutgeneratingfunctionology.igp import *
from minimalfunctioncache.system import sys_info
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
# sys_info is an object that contains knows how to read data from the minimal function cache.

# current defaults
max_bkpts = sys_info.max_bkpts
tol = sys_info.global_default_tol
solver_mode_default = "full"
row_processing_order_default = "lex"
cut_scoring_method_default = "parallelism"

class UnsetData(Exception):
    pass

class ParamaterizationError(Exception):
    pass

class abstractCutScore:
    r"""
    abstract class for cut optimization objective functions, aka cut scoring metrics. 
    """
    @classmethod
    def __init__(cls, **kwrds):
        pass

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

class cutScore:
    """
    cutScore is objective function used in the cutOptimzationProblem.

    cutScore is fixed once data from the problem is fixed.
    """
    @staticmethod
    def __classcall__(cls, name=None, **kwrds):
        r"""
        input normalization of class
        """
        if name == "parallelism" or name is None:
            return super().__classcall__(cls, cut_score=Parallelism)
        if name == "steepest_direction" or name is None:
            return super().__classcall__(cls, cut_score=SteepestDirection)
        if issubclass(name, abstractCutScore):
            return super().__classcall__(cls, cut_score=name, **kwrds)
        else:
            raise TypeError("Use a predefined cut scoring method or use custom instance of abstractCutScore.")

    def __init__(self, cut_score, **kwrds):
        r"""
        Initialize the paramatrized cut scoring function. 
        
        Maximization or minimlzation should be specified upon initlization.
        
        Data used with cutScore is managed by cutGenerationSolverBase.
        """
        self._cut_score = cut_score(**kwrds)
        super().__init__()
        self._MIP_objective = None
        self._MIP_row = None
        self._sage_to_solver_type = None
        if "obj_type" in kwrds.keys():
            self._cut_obj_type = kwrds.keys()["obj_type"]
        else:
            self._cut_obj_type = "max"
    
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

    def set_MIP_objective(self, new_objective):
        """
        Use reduced costs of the basis relaxation here. 

        This should be a vector of length n-m for the problem. 
        """
        self._MIP_objective = new_objective

    def get_MIP_objective(self):
        return self._MIP_objective

    def get_MIP_row(self):
        """
        for fixed i,bar_i - (bar a_i)^T x_N 
        """
        return self._MIP_row

    def set_MIP_row(self, new_row):
        self._MIP_row = new_row

    def _get_sage_to_solver_type(self):
        """
        Defined in the cutGenerationSolver. This method is only indeded to be set and unset by solving
        routines in cutGenerationSolver. 
        """
        return self._sage_to_solver_type

    def _set_sage_to_solver_type(self, new_conversion):
        self.__sage_to_solver_type = new_conversion

    def cut_obj_type(self):
         self._cut_obj_type

class Parallelism(abstractCutScore):
    """
    Normalized cut parallelism score.
    """
    def cut_score(cls, cut, field=None):
        obj_norm = vector(cls._MIP_objective).norm()
        cut_norm = vector(cut).norm()
        dot_product = vector(cls._MIP_objective).row()*vector(cut).column()
        return dot_product[0]/(obj_norm*cut_norm)
        
class SteepestDirection(abstractCutScore):
    """
    Steepest direction score. 
    """
    def cut_score(cls, cut, field=None):
        dot_product = vector(cls._MIP_objective).row()*vector(cut).column()
        return dot_product[0]
        

class cutGenerationDomain:
    r"""PiMin<=k, possibly with paramaterized constraints.
    """
    def __init__(self, k):
        self._num_bkpts = k
        self._PiMin = None
        self._constraints = []

    def __repr__(self):
        return "The Space of Con. Minimal Functions with at most {} breakpoints and constrained by {}.".format(self._num_bkpts, self._constraints)
    # def is_cell_description(self):
        # pass
    # def is_manifold_description(self):
        # pass
    
    def _load_PiMin(self):
        if k > max_bkpts: # and len(self.possible_bkpt_params) == 0: 
           raise ValueError("The minimal function cache requested has not been computed. Please compute the minimal funciton cache for {}.".format(k))
    
    def add_constraint(self, con):
        self._constraints.append(con)
        
    def get_cells(self):
        r"""iterate over all cells (either explicilty defined or infered) with intersected constraints. Outputs a BSA.
        """
        if self._PiMin is None:
            self._load_PiMin()
        cells = self._PiMin.get_semialgebraic_sets()
        for cell in cells:
            yield cell.intersection(self._constraints)

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
        return value_nnc_polyhedron(bkpt, f_index) # write a version which includes makes a consistent number of parameters. 


class cutGenerationSolverBase:
    r"""
    A base class for interfacing with solvers. cutGenerationDomains can be expressed a finite collection of semialgebraic cells. 
    This has methods translate semialgebraic cell or manifold descriptions of data from the cutGenerationDomain to some external solver.
    The solve and __init__ methods are written to be compatiable with the cut generation problem  and provide an optimal solution to the problem 
    provided that the solver and data conversion methods have been written. 
    
    cutGenerationSovler infers the correct use of cutGenerationDomain based on provided options from the cutGenerationProblem. 
    
    The cutGenerationProblem options are listed below.
    
    Option: row selection: all_rows mathcal I = {i \in B : overline{b_i} \not \in ZZ}, lex_row = min  mathcal I, random row, rand(mathcal I), subset (any subset of mathcal I)
    Option: algorithm = full space; bkpt as param: full (all combinatorial data, if needed k<n-m),  max, lex (selects lexicographically first parameter data k< n-m), or rand (randomly select parameter data)
    Option: num_bkpt = full space: k <= max_bkpts; bkpt as param: full, lex rand, k <= n-m; max, k = n-m
    Option: cut_score = parallelism, steepestdirection, scip, or custom
    Option: multithread = notImplemented
    """
    def __init__(self, algorithm=None, cut_score=None, num_bkpt=None, row_selection=None, multithread=False):
        if algorithm is None:
            self. _algorithm = "bkpt as param, full"
        self._cut_score = cut_score
        self._cut_score.set_sage_to_solver_type(self.sage_to_solver_type)
        self._algorithm = algorithm
        self._num_bkpt = num_bkpt
        

    def solve(self, MIP, **paramaterized_problem_solver_options):
        r"""Solves the paramaterized problem. 
        
        Interperts the options and calls the correct solving algorithm. 
        
        Passes any instructions to the underlying solver.
        
        """
        # for now assuming problem is in the form of min
        
        if self._algorithm = "full":
            result = self._algorithm_full_space(MIP)
        elif self._algorithm =  "bkpt as param, full":
            result = self._algorithm_bkpt_as_param_full(MIP)
        return result

    def _algorithm_full_space(self, MIP):
        r"""
        """
        pass
        
        # for subdomain in self._cut_gen_domain.get_full_cells():
            # objective_fun = self._cut_score # the call method here is a function from R^{2n} (parameter space) to R
            # subdomain_solver_constraints = self.write_nonlinear_constraints_from_bsa_for_solver(subdomain)
            # subdomain_problem_val, subdomain_problem_solution = self.solver_nonlinear_solve(subdomain_solver_constraints, objective_fun, min_or_max, **solver_options)                

    def _algorithm_bkpt_as_param_full(self, MIP):
        r"""
        """
        pass

    @staticmethod
    def write_linear_constraints_from_bsa_for_solver(self, bsa): # think about aspects of exactness; 
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
    

class sciPYCutGenSolvers(cutGenerationSolverBase):
    @staticmethod
    def write_linear_constraints_from_bsa_for_solver(self, bsa, epsilon=10**-9): # think about aspects of exactness; 
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """ 
        pass

    @staticmethod       
    def write_nonlinear_constraints_from_bsa(self, bsa, epsilon=10**-9):
        r"""
        Given a BSA with non linear constraints, converts the bsa object into a format that the underlying solver can use.
        
        Treats p(x) < c as p(x) + epsilon <= c for all epsilon>0.
        """
        nonlinear_constraints = []
        for con in bsa.eq_poly():
            def con_conv_fun(array_like):
                input_map = {con.parent().gens()[i]: array_like[i] for i in range(con.parent().ngens()}
                return np.array([con.subs(input_map)])
            nonlinear_constraints.append(NonlinearConstraint(con_conv_fun, 0,0))
        for con in bsa.le_poly():
            def con_conv_fun(array_like):
                input_map = {con.parent().gens()[i]: array_like[i] for i in range(con.parent().ngens()}
                return np.array([con.subs(input_map)])
            nonlinear_constraints.append(NonlinearConstraint(con_conv_fun, -np.inf, 0))
        for con in bsa.lt_poly():
            def con_conv_fun(array_like):
                input_map = {con.parent().gens()[i]: array_like[i] for i in range(con.parent().ngens()}
                return np.array([con.subs(input_map)+epsilon])
            nonlinear_constraints.append(NonlinearConstraint(con_conv_fun, -np.inf, 0))        
        # what needs to happen is composion of maps. we have map that that goes from cords lambda1,....,lambdak, gamma1,....,gammak
        return nonlinear_constraints


    @staticmethod
    def solver_linear_solve(self, constraints, objective, min_or_max,  **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b.      
        
        Should return optimal objective value, optimal objective solution.
        """
        raise NotImplementedError

    @staticmehtod
    def solver_nonlinear_solve(self, constraints, objective, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polynomial and at least 1 p_i has degree larger than 1.  
        """
        minimize(objective, constraints) # think about this...


    @staticmethod
    def sage_to_solver_type(self, sage_field_element, field=None):
        return float(sage_field_element)


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
