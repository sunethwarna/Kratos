from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
# from analysis_stage import AnalysisStage

# Importing derived classes
from potential_flow_analysis import PotentialFlowAnalysis

# Import pycompss
# from pycompss.api.task import task
# from pycompss.api.api import compss_wait_on
# from pycompss.api.parameter import *

# Import exaqute
# from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
# get_value_from_remote is the equivalent of compss_wait_on
# in the future, when everything is integrated with the it4i team, putting exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC

# Import variables class
from cmlmc_utilities import StatisticalVariable

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

class MonteCarloAnalysis(PotentialFlowAnalysis):
    """Main analyis stage for Monte Carlo simulations"""
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(input_model,input_parameters)
        # self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _GetSimulationName(self):
        return "Monte Carlo Analysis"

    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        Mach = self.sample[0]
        a_infinity = 340 # [m/s] velocity of sound at infinity
        alpha =  self.sample[1]
        v_norm = Mach * a_infinity
        velocity = [v_norm*np.cos(alpha),v_norm*np.sin(alpha),0]
        boundary_processes = self.project_parameters["processes"]["boundary_conditions_process_list"]
        problem_name=self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        for i in range(0,boundary_processes.size()):
            python_module = boundary_processes[i]["python_module"].GetString()
            if python_module == "apply_far_field_process":
                self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["velocity_infinity"].SetVector(velocity)
        # self.project_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(problem_name+'_M'+str(Mach)+'_A'+str(alpha))


##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################


'''
function generating the random sample
'''
def GenerateSample():
    sample = []
    mean_Mach = 0.3
    std_deviation_Mach = 0.01
    number_samples = 1
    sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
    mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
    std_deviation_angle_attack = np.deg2rad(0.1)
    sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
    if sample[0] >= 1.0 or sample[0] <= 0.0 :
        raise Exception ("stochastic Mach number computed > 1 or < 0")
    return sample


'''
function evaluating the QoI of the problem
'''
def EvaluateQuantityOfInterest(simulation):
    """here we evaluate the QoI of the problem: the lift coefficient"""
    Q = simulation._GetSolver().main_model_part.GetValue(KratosMultiphysics.FRICTION_COEFFICIENT)
    return Q

'''
function executing the problem
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        QoI         : Quantity of Interest
'''
@ExaquteTask(returns=1)
def ExecuteMonteCarlo_Task(pickled_model, pickled_parameters):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    sample = GenerateSample()
    simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    return QoI


'''
function serializing the model and the parameters of the problem
the idea is the following:
first from Model/Parameters Kratos object to StreamSerializer Kratos object
second from StreamSerializer Kratos object to pickle string
third from pickle string to StreamSerializer Kratos object
fourth from StreamSerializer Kratos object to Model/Parameters Kratos object
input:
        parameter_file_name   : path of the Project Parameters file
output:
        pickled_model      : model serialized
        pickled_parameters : project parameters serialized
'''
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters_Task(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    fake_sample = GenerateSample() # fake sample generated just to do the serialization
    simulation = MonteCarloAnalysis(model,local_parameters,fake_sample)
    simulation.Initialize()
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2)
    return pickled_model, pickled_parameters


if __name__ == '__main__':

    parameter_file_name = "../tests/CompressiblePotentialFlowTest/parameters_potential_naca_mesh0.json"

    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    print("\n############## Serialization completed ##############\n")

    '''estimate the number of samples'''
    number_samples = 3

    QoI = StatisticalVariable(0) # number of levels = 0 (we only have one level), needed using this class
    '''to exploit StatisticalVariable UpdateOnePassMeanVariance function we need to initialize a level 0 in values, mean, sample variance and second moment
    and store in this level the informations'''
    QoI.values = [[] for _ in range (1)]
    QoI.mean = [[] for _ in range (1)]
    QoI.second_moment = [[] for _ in range (1)]
    QoI.sample_variance = [[] for _ in range (1)]

    for instance in range (0,number_samples):
        QoI.values[0].append(ExecuteMonteCarlo_Task(pickled_model,pickled_parameters))

    '''Compute mean, second moment and sample variance'''
    for i_sample in range (0,number_samples):
        QoI.UpdateOnepassMeanVariance(0,i_sample)

    QoI = get_value_from_remote(QoI)
    print("\nMonte Carlo mean estimator lift coefficient = ",QoI.mean)

