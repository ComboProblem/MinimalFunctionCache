import tomllib
import logging
from pyscipopt import Model
from cutGen.optimal_cut_generation import OptimalCut


test_model = Model()

def scip_parameter_parser(model, path_to_scip_config_toml):
    """configures the model"""
    pass
    
def SciPy_paramater_parser(path_to_SciPy_config_toml):
    pass
    
def cutGenProb_parameter_parser(path_to_cutGenProb_config_toml, configured_cutGenProblemSolver):
    """
    Input - 
    Output - configured optimal cut
    """
    pass
    


# parse experimental setting

# create jobs by sending job to runner