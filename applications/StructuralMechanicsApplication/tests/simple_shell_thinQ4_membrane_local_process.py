from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):

        self.model_part = model_part[params["model_part_name"].GetString()]
        self.params = params
        
    def ExecuteInitialize(self):
        pass
        
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        
        #determine max y disp
        dy_max = 0
        for node in self.model_part.Nodes:
            dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
            if(dy > dy_max):
                dy_max = dy
        self.assertAlmostEqual(dy_max,0.05, places = 6)
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
