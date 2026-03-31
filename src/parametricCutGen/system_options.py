class sys_info:
    r"""
    Class for mutable solver options. Normalizes inputs for experiements.  
    """
    def __init__(self):
        self.max_bkpts = None
        self.use_k_bkpts = None
        self.global_default_tol = None
        self.solver_mode = None
        self.cut_scoring_method = None
        self.sys_info.cut_problem_solver = None