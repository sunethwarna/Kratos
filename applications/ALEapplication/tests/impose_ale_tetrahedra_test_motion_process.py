from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *


def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeALETetrahedraTestMotionProcess(Model, settings["Parameters"])

class ImposeALETetrahedraTestMotionProcess(Process):
    def __init__(self, Model, settings):
        Process.__init__(self)
        
        self.model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        for node in self.model_part.Nodes:
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)
        
    def ExecuteInitializeSolutionStep(self):
        time = self.model_part.ProcessInfo[TIME]
        xc = 0.375
        yc = 0.250
        pi = 3.141592653589793
        ut = 0.5 * sin(pi * time)
        phi = 0.5 * sin(3.0 * pi * time)
        c = cos(phi)
        s = sin(phi)
        for node in self.model_part.Nodes:
            uz = 1.0 * time
            rx = node.X0 - xc
            ry = node.Y0 - yc
            ur = xc + c * rx - s * ry - node.X0
            vr = yc + s * rx + c * ry - node.Y0
            node.SetSolutionStepValue(MESH_DISPLACEMENT_X, ut + ur)
            node.SetSolutionStepValue(MESH_DISPLACEMENT_Y, vr)
            node.SetSolutionStepValue(MESH_DISPLACEMENT_Z, uz)
