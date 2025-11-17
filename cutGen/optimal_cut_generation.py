from cutgeneratingfunctionology.igp import *
from minimalfunctioncache.system import sys_info
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from pyscipopt import Model, Sepa
# sys_info is an object that contains knows how to read data from the minimal function cache.

# current defaults
max_bkpts = sys_info.max_bkpts
use_k_bkpts = sys_info.use_k_bkpts
tol = sys_info.global_default_tol
solver_mode = sys_info.solver_mode
algo  =  sys_info.algorithm
sys_threading = sys_info.multithreading
cut_scoring_method = sys_info.cut_scoring_method
cut_problem_solver = sys_info.cut_problem_solver

def find_f_index(min_pwl):
    min_pwl.end_points().index(find_f(min_pwl))

class UnsetData(Exception):
    pass

class abstractCutScore:
    r"""
    Abstract class for cut optimization objective functions aka cut scores.
    Named after huersticis used to evalaute a cuts effectiveness.
    """
    @classmethod
    def __init__(cls, **kwrds):
        pass

    def cut_score(cls, cut, mip_obj):
        r"""
        A(n) (assumed to be) smooth function from cutSpace to RR where cutSpace = {(\pi(bar a_ij))_{j\in N} : pi in PiMin}. 

        Suppose that R is a ring such that either QQ subseteq R subseteq RR or elements of R 
        coherese to a ring R' wiht QQ subseteq R' subseteq RR.

        cut_score should use sagemath types to ensure generating a seperating cut.
        
        input: cut: A list of elements of R
        mip_obj: A list of elements of R

        output: an element of R

        EXAMPLES:
        """
        raise NotImplementedError

    def cut_score_grad(cls, cut, mip_obj):
        r"""
        """
        raise NotImplementedError
    
    def cut_score_hess(cls, cut, mip_obj):
        r"""
        """
        raise NotImplementedError

class cutScore:
    """
    cutScore is objective function used in the cutOptimzationProblem.
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

    def __init__(self, cutscore, **kwrds):
        r"""
        Initialize the paramatrized cut scoring function. 
        
        Maximization or minimlzation should be specified upon initlization.
        
        Data used with cutScore is managed by cutGenerationSolverBase.
        """
        self._cut_score = cutscore(**kwrds)
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
        pi = piecewise_function_from_breakpoints_and_values(bkpt + [1], val + [0])
        row_data = self.get_MIP_row()
        sage_cut = [pi(fractional(QQ(bar_a_ij))) for bar_a_ij in row_data]
        sage_mip_obj =  [QQ(bar_cj) for bar_cj in self._MIP_objective]
        result = self.get_sage_to_solver_type()(self._cut_score.cut_score(sage_cut, sage_mip_obj))
        return result

    def set_MIP_obj(self, new_objective):
        """
        Use reduced costs of the basis relaxation here. 

        This should be a vector of length n-m for the problem. 
        """
        self._MIP_objective = new_objective

    def get_MIP_obj(self):
        return self._MIP_objective

    def get_MIP_row(self):
        """
        for fixed i,bar_i - (bar a_i)^T x_N 
        """
        return self._MIP_row

    def set_MIP_row(self, new_row):
        self._MIP_row = new_row

    def get_sage_to_solver_type(self):
        """
        Defined in the cutGenerationSolver. This method is only indeded to be set and unset by solving
        routines in cutGenerationSolver. 
        """
        return self._sage_to_solver_type

    def set_sage_to_solver_type(self, new_conversion):
        self._sage_to_solver_type = new_conversion

    def cut_obj_type(self):
         self._cut_obj_type

class Parallelism(abstractCutScore):
    """
    Normalized cut parallelism score.
    """
    def cut_score(cls, cut, mip_obj):
        obj_norm = vector(mip_obj).norm()
        cut_norm = vector(cut).norm()
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return (dot_product[0]/(obj_norm*cut_norm))[0]
        
class SteepestDirection(abstractCutScore):
    """
    Steepest direction score. 
    """
    def cut_score(cls, cut, mip_obj):
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return dot_product[0]


class cutGenerationProblem:
    r""" 
    A base class for interfacing with solvers and solvoing a row cut generation problem.
    
    
    cutGenerationSovler infers the correct use of cutGenerationDomain based on provided options from the cutGenerationProblem. 
                
    The cutGenerationProblem options are listed below.
    
    Option: row selection: all_rows mathcal I = {i \in B : overline{b_i} \not \in ZZ}, lex_row = min  mathcal I, random row, rand(mathcal I), subset (any subset of mathcal I)
    Option: row algorithm = full space; bkpt as param: full (all combinatorial data, if needed k<n-m),  max, lex (selects lexicographically first parameter data k< n-m), or rand (randomly select parameter data)
    Option: num_bkpt = full space: k <= max_bkpts; bkpt as param: full, lex rand, k <= n-m; max, k = n-m
    Option: cut_score = parallelism, steepestdirection, scip, or custom
    Option: multithread = notImplemented
    """
    def __init__(self, algorithm=None, cut_score=None, num_bkpt=None, solver=None, multithread=False):
        if algorithm is None:
            self. _algorithm = "bkpt as param, full"
        self._cut_score = cut_score
        self._solver = solver
        self._cut_score.set_sage_to_solver_type(self._solver.sage_to_solver_type)
        self._algorithm = algorithm
        self._num_bkpt = num_bkpt
        self._cut_space = None
        

    def solve(self, binvarow, binvc, f):
        r"""Solves the paramaterized problem. 
        
        Interperts the options and calls the correct solving algorithm. 
        
        Passes any instructions to the underlying solver.
        
        """
        # assume MIP is a scip model; really we should be passing in and LP relaxation with variable infomation here.
        # The cut generation problem 
        if self._algorithm == "full":
            result = self._algorithm_full_space(binvarow, binvc, f)
        elif self._algorithm == "bkpt as param, full":
            result = self._algorithm_bkpt_as_param_full(binvarow, binvc, f)
        return result

    def _algorithm_full_space(self, binvarow, binvc, f):
        r"""
        Solves the problem given a row of B^-1A and the reduced costs 
        """
        self._cut_score.set_MIP_row(binvarow)
        self._cut_score.set_MIP_obj(binvc)
        def objective_fun(params):
            return self._cut_score(params) #lambda x: self._cut_score(x)
        # if max_or_min == "max":
        best_value = -1*np.inf
        # else: #To do implement minimze options
        #     raise NotImplementedError
        result =  None
        if self._cut_space is None: # load the semi algebraic descriptions. 
             if self._num_bkpt > max_bkpts:
                raise ValueError("The Minimal Functions Cache for {} breakpoints requested has not been computed.".format(self._num_bkpt))
             self._cut_space = PiMinContContainer(self._num_bkpt)
        for b, v in self._cut_space.get_rep_elems():
            # f is a bkpt when pi has a finite number of bkpts.
            pi_test = piecewise_function_from_breakpoints_and_values(b+[1], v+[0])
            bsa_f_index = find_f_index(pi_test)
            subdomain_with_f_constraint = bsa_of_rep_element(b, v)
            lambda_f_index = subdomain_with_f_constraint.polynomial_map()[0].parent().gens()[bsa_f_index]
            lhs =  lambda_i - QQ(f)
            x0 = b+v # pick inital guess to be a point known in the cell
            subdomain_with_f_constraint.add_polynomial_constraint(lhs, operator.eq)
            # TODO: Use bkpt to find if 
            subdomain_solver_constraints = self._solver.write_nonlinear_constraints_from_bsa(subdomain_with_f_constraint)
            subproblem_result = self._solver.nonlinear_solve(objective_fun, x0, subdomain_solver_constraints)   
            if subproblem_result[0] is True:
                if best_value < subproblem_result[1]:
                    best_value = subproblem_result[1]
                    result = subproblem_result
        # if result is None, the solver has failed to find any meaninful result. There should always be a result 
        # since gmic(f) bounds the maximization problem from above.
        # TODO: Add something to reterive full status of subproblem solution.
        # result [0
        bkpt_result = result[2][self._num_bkpt/2:]
        vals_result = result[2][:self._num_bkpt/2]
        pi_p = piecewise_function_from_breakpoints_and_values(bkpt_result+[1],vals_result+[0])
        return pi_p

    def _algorithm_bkpt_as_param_full(self, rbinvarow, binvc, f):
        pass

    def _algorithm_custom(self, relaxed_row, mip_obj, f):
        raise NotImplementedError

class abstractCutGenProblemSolverInterface:
    r"""
    Interfaces types from ``cutgeratingfunctionolgy`` to a specified solver.
    """
    def __init__():
        pass
    @staticmethod
    def write_linear_constraints_from_bsa(bsa):
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """     
        raise NotImplementedError

    @staticmethod
    def write_nonlinear_constraints_from_bsa(bsa):
        r"""
        Given a BSA with non linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError


    @staticmethod
    def linear_solve(constraints, objective,  **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b. 
        
        Should return optimal objective value, optimal objective solution, solver success, and solver_output
        """
        raise NotImplementedError
    
    @staticmethod
    def nonlinear_solve(constraints, objective, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polynomial and at least 1 p_i has degree larger than 1.  
        """
        raise NotImplementedError


    @staticmethod
    def sage_to_solver_type(sage_ring_element):
        r"""
        
        """
        raise NotImplementedError
    

class scipyCutGenProbelmSolverInterface(abstractCutGenProblemSolverInterface):
    """
    Interfaces types and objects from ``cutgeneratingfunctiology`` to scipy. 
    """
    @staticmethod
    def write_linear_constraints_from_bsa(bsa, epsilon=10**-9): # think about aspects of exactness; 
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """ 
        pass

    @staticmethod       
    def write_nonlinear_constraints_from_bsa(bsa, epsilon=10**-9):
        r"""
        Given a BSA with nonlinear constraints, converts into an equivlent set of nonlinear constraints for scipy. 
        
        Treats p(x) < c as p(x) + epsilon <= c for all epsilon>0.
        """
        # Read polynomial parent to define map from BSA ambident to space to solver ambient space.
        # input_map is equivlent to the identiy map. 
        
        # TODO: since BSA is only has polynomial constraints, explicitly define the gradient and hessian 
        # of each constraint. 
        # Either do some type of ring theoretic dual number extension,
        # Add an auto diff feature or 
        # or write scripts to explicitly computer the derivatives since polynomial derivatives are not that hard.
        nonlinear_constraints = []
        for polynomial in bsa.eq_poly():
            def poly(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)])
            nonlinear_constraints.append(NonlinearConstraint(poly, 0,0))
        for polynomial in bsa.le_poly():
            def poly(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)])
            nonlinear_constraints.append(NonlinearConstraint(poly, -np.inf, 0))
        for polynomial in bsa.lt_poly():
            def poly(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)+epsilon])
            nonlinear_constraints.append(NonlinearConstraint(poly, -np.inf, 0))        
        return nonlinear_constraints


    @staticmethod
    def lp_solve(objective, constraints, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b.      
        
        Should return optimal objective value, optimal objective solution.
        """
        raise NotImplementedError

    @staticmethod
    def nonlinear_solve(objective, x0, cons, **solver_options):
        r"""
        Given converted constraints and an objective function that is compitable with the solver, 
        use scipy minimize to ...
        """
        
        result =  minimize(objective, x0, constraints=cons) # think about this...
        return result.success, result.fun, result.x, result

    @staticmethod
    def sage_to_solver_type(sage_ring_element):
        """
        scipy supports inputs of floats, convert sage field element to its equivlant numerical (python) floating point value.
        """
        
        # be lazy and assume the_ring_element is something the converts to a rational number (or can be put into a floating point approximation).
        # this is a point where we lose the exactness of sage. 
        return float(sage_ring_element)




# Following SCIP GMIC tutorial; naturally since this is a genearlization of this type of seperator business. 
# This seems the easist way to interface wiht SCIP. I will use the problem as an interface dispatching
# options and processing datta into optimal cuts. In some shense, the problem needs
# pramaterize data, then return cuts. 
# The solver returns paramaterized data on from row MIP data
# so the problem is responsible for putting MIP row data into a sover compitable form.
# problem should serve as the interface between MIP solver and the cut Problme Solver. 
# theCutProblem needs to take MIP data and put in into sage exact types. 
# Baiscly, we want the paramaterization maps to be exact. So after solver and we've produced
# a cut that is indeed and intersection cut for the numeridcal data presented. 
# The therotical gareuntee of producing intersection cuts is what we are after. And we want from the prespective for the MIP's
# solvers numerics to be the only worry for solving the MIP as the cut produced from the cut problem is exact (from a rational prespective)
# once the cut is converted to the MIP, depending on how the floating point arithmatic happens, we might lose this exact guantree, but
# at least we will attempt to take steps to ensure that property is perserved. 


# Modifying the class from SCIP tutorial. All that is needed here is to rewrite the cutelem defintion to be pi_p. I can insert 
# any GEMIC funciton


class OptimalCut(Sepa):
    def __init__(self):
        self.ncuts = 0
        self.cgp = cutGenerationProblem(algorithm=algo, cut_score=cut_scoring_method, num_bkpt=use_k_bkpts, multithread=sys_threading, solver=cut_problem_solver)
    def getOptimalCutFromRow(self, cols, rows, binvrow, binvarow, primsol, pi_p):
        """ Given the row (binvarow, binvrow) of the tableau, computes optimized cut

        :param primsol:  is the rhs of the tableau row.
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A
        :param 1dPWL:    a minimal funciton element of PiMin<=k with pi_p(primsol)=1

        The interesection cut is given by
         sum(pi_p(a_j) x_j, j in J_I) \geq 1
        where J_I are the integer non-basic variables and J_C are the continuous.
        f_0 is the fractional part of primsol
        a_j is the j-th coefficient of the row and f_j its fractional part
        Note: we create -% <= -f_0 !!
        Note: this formula is valid for a problem of the form Ax = b, x>= 0. Since we do not have
        such problem structure in general, we have to (implicitly) transform whatever we are given
        to that form. Specifically, non-basic variables at their lower bound are shifted so that the lower
        bound is 0 and non-basic at their upper bound are complemented.
        """

        # initialize
        cutcoefs = [0] * len(cols)
        cutrhs = 0

        # get scip
        scip = self.model


        # Generate cut coefficients for the original variables
        for c in range(len(cols)):
            col = cols[c]
            assert col is not None
            status = col.getBasisStatus()

            # Get simplex tableau coefficient
            if status == "lower":
                # Take coefficient if nonbasic at lower bound
                rowelem = binvarow[c]
            elif status == "upper":
                # Flip coefficient if nonbasic at upper bound: x --> u - x
                rowelem = -binvarow[c]
            else:
                # variable is nonbasic free at zero -> cut coefficient is zero, skip OR
                # variable is basic, skip
                assert status == "zero" or status == "basic"
                continue

            # Integer variables
            if col.isIntegral():
                # warning: because of numerics cutelem < 0 is possible (though the fractional part is, mathematically, always positive)
                # However, when cutelem < 0 it is also very close to 0, enough that isZero(cutelem) is true, so we ignore
                # the coefficient (see below)
                cutelem = float(pi_p(fractional(QQ(rowelem)))) #keep types correct
            else:
                # Continuous variables
                # if rowelem < 0.0:
                    # -sum(a_j*f_0/(1-f_0) x_j      , j in J_C s.t. a_j  <   0) >= f_0.
                    # cutelem = rowelem * ratiof0compl
                # else:
                    #  sum(a_j x_j,                 , j in J_C s.t. a_j >=   0) -
                    # cutelem = -rowelem
                pass # IPs don't have cont. variables, which is what I'm writing for right  now. 
            # cut is define when variables are in [0, infty). Translate to general bounds
            if not scip.isZero(cutelem):
                if col.getBasisStatus() == "upper":
                    cutelem = -cutelem
                    cutrhs += cutelem * col.getUb()
                else:
                    cutrhs += cutelem * col.getLb()
                # Add coefficient to cut in dense form
                cutcoefs[col.getLPPos()] = cutelem

        # Generate cut coefficients for the slack variables; skip basic ones
        for c in range(len(rows)):
            row = rows[c]
            assert row != None
            status = row.getBasisStatus()

            # free slack variable shouldn't appear
            assert status != "zero"

            # Get simplex tableau coefficient
            if status == "lower":
                # Take coefficient if nonbasic at lower bound
                rowelem = binvrow[row.getLPPos()]
                # But if this is a >= or ranged constraint at the lower bound, we have to flip the row element
                if not scip.isInfinity(-row.getLhs()):
                    rowelem = -rowelem
            elif status == "upper":
                # Take element if nonbasic at upper bound - see notes at beginning of file: only nonpositive slack variables
                # can be nonbasic at upper, therefore they should be flipped twice and we can take the element directly.
                rowelem = binvrow[row.getLPPos()]
            else:
                assert status == "basic"
                continue

            # if row is integral we can strengthen the cut coefficient
            if row.isIntegral() and not row.isModifiable():
                # warning: because of numerics cutelem < 0 is possible (though the fractional part is, mathematically, always positive)
                # However, when cutelem < 0 it is also very close to 0, enough that isZero(cutelem) is true (see later)
                cutelem = float(pi_p(fractional(QQ(rowelem))))
            else:
                # Continuous variables
                # if rowelem < 0.0:
                    # -sum(a_j*f_0/(1-f_0) x_j      , j in J_C s.t. a_j  <   0) >= f_0.
                #     cutelem = rowelem * ratiof0compl
                # else:
                    #  sum(a_j x_j,                 , j in J_C s.t. a_j >=   0) -
                    # cutelem = -rowelem
                pass

            # cut is define in original variables, so we replace slack by its definition
            if not scip.isZero(cutelem):
                # get lhs/rhs
                rlhs = row.getLhs()
                rrhs = row.getRhs()
                assert scip.isLE(rlhs, rrhs)
                assert not scip.isInfinity(rlhs) or not scip.isInfinity(rrhs)

                # If the slack variable is fixed, we can ignore this cut coefficient
                if scip.isFeasZero(rrhs - rlhs):
                  continue

                # Unflip slack variable and adjust rhs if necessary: row at lower means the slack variable is at its upper bound.
                # Since SCIP adds +1 slacks, this can only happen when constraints have a finite lhs
                if row.getBasisStatus() == "lower":
                    assert not scip.isInfinity(-rlhs)
                    cutelem = -cutelem

                rowcols = row.getCols()
                rowvals = row.getVals()

                assert len(rowcols) == len(rowvals)

                # Eliminate slack variable: rowcols is sorted: [columns in LP, columns not in LP]
                for i in range(row.getNLPNonz()):
                    cutcoefs[rowcols[i].getLPPos()] -= cutelem * rowvals[i]

                act = scip.getRowLPActivity(row)
                rhsslack = rrhs - act
                if scip.isFeasZero(rhsslack):
                    assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                    cutrhs -= cutelem * (rrhs - row.getConstant())
                else:
                    assert scip.isFeasZero(act - rlhs)
                    cutrhs -= cutelem * (rlhs - row.getConstant())

        return cutcoefs, cutrhs

   def sepaexeclp(self):
        result = SCIP_RESULT.DIDNOTRUN
        scip = self.model

        if not scip.isLPSolBasic():
            return {"result": result}
    
        # get LP data
        cols = scip.getLPColsData()
        rows = scip.getLPRowsData()

        # exit if LP is trivial
        if len(cols) == 0 or len(rows) == 0:
            return {"result": result}

        result = SCIP_RESULT.DIDNOTFIND

        # get basis indices
        basisind = scip.getLPBasisInd()

        # For all basic columns (not slacks) belonging to integer variables, try to generate a gomory cut
        for i in range(len(rows)):
            tryrow = False
            c = basisind[i]

            if c >= 0:
                assert c < len(cols)
                var = cols[c].getVar()

                if var.vtype() != "CONTINUOUS":
                    primsol = cols[c].getPrimsol()
                    assert scip.getSolVal(None, var) == primsol

                    if 0.005 <= scip.frac(primsol) <= 1 - 0.005:
                        tryrow = True

            # generate the cut!
            if tryrow:
                # get the row of B^-1 for this basic integer variable with fractional solution value
                binvrow = scip.getLPBInvRow(i)

                # get the tableau row for this basic integer variable with fractional solution value
                binvarow = scip.getLPBInvARow(i)
               
                # get current reduced costs for objective evaluation.
                costs = [scip.getColRedCost(j) for j in range(len(cols)) if j not inbasisind]

                cgf = self.cgp.solve(binvarow, costs, primsol) # produce an optimal cgf
                
                cutcoefs, cutrhs = self.getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)

                # add cut
                cut = scip.createEmptyRowSepa(self, "gmi%d_x%d"%(self.ncuts,c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
                scip.cacheRowExtensions(cut)

                for j in range(len(cutcoefs)):
                    if scip.isZero(cutcoefs[j]): # maybe here we need isFeasZero
                        continue
                    scip.addVarToRow(cut, cols[j].getVar(), cutcoefs[j])

                if cut.getNNonz() == 0:
                    assert scip.isFeasNegative(cutrhs)
                    return {"result": SCIP_RESULT.CUTOFF}


                # Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound changes)
                # the latter cuts will be handled internally in sepastore.
                if cut.getNNonz() == 1 or scip.isCutEfficacious(cut):

                    # flush all changes before adding the cut
                    scip.flushRowExtensions(cut)

                    infeasible = scip.addCut(cut, forcecut=True)
                    self.ncuts += 1

                    if infeasible:
                       result = SCIP_RESULT.CUTOFF
                    else:
                       result = SCIP_RESULT.SEPARATED
                scip.releaseRow(cut)

        return {"result": result}
