import tomllib
import logging
import os
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut

experiments_logger = logging.getLogger(__name__)

def scip_parameter_parser(model, path_to_scip_config_toml):
    """configures the model"""
    pass
    
def SciPy_paramater_parser(path_to_SciPy_config_toml):
    pass
    


def experiment(path_to_scip_config_toml,path_to_SciPy_config_toml, experiment_parameters):
    experimental_model = Model()
    cutGen = 
    scip_parameter_parser( ,path_to_scip_config_toml)
    seap = OptimalCut(algorithm_name=experiment_parameters["algorithm"], cut_score=experiment_parameters["cutScore"], num_bkpt=experiment_parameters["num_bkpt"], multithread=False, prove_seperator=experiment_parameters["prove_seperator"], show_proof = False, epsilon=10**-7, M = 10**7)

class Experiement:
    """
    Parametrically defined experiments. 
    
    Parameters: Algorithm
    If ne
    """
    
    def __init__(**parameters):
        """
        Checks the experiemental 
        """
        self._experiment_model = Model()
    def run_experiment():
        self._experiment_model
        
    def write_results():
        pass

    def get_logs():
        pass




