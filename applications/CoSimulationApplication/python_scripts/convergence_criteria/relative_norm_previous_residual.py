from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_criteria import CoSimulationConvergenceCriteria

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import numpy as np
from numpy import linalg as la

def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return RelativeNormPreviousResidualConvergenceCriteria(settings)

class RelativeNormPreviousResidualConvergenceCriteria(CoSimulationConvergenceCriteria):
    def __init__(self, settings):
        super(RelativeNormPreviousResidualConvergenceCriteria, self).__init__(settings)

        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()
        self.label = self.settings["label"].GetString()

    def IsConverged(self, residual, current_data):
        res_norm = la.norm(residual)
        norm_new_data = la.norm(current_data)

        if norm_new_data < 1e-15:
            norm_new_data = 1.0 # to avoid division by zero

        abs_norm = res_norm / np.sqrt(residual.size)
        rel_norm = res_norm / norm_new_data

        is_converged = abs_norm < self.abs_tolerance or rel_norm < self.rel_tolerance

        if self.echo_level > 1:
            info_msg  = 'Convergence '

            if self.label != "":
                info_msg += 'for "{}": '.format(self.label)

            if is_converged:
                info_msg += colors.green("ACHIEVED")
            else:
                info_msg += colors.red("NOT ACHIEVED")
            cs_tools.cs_print_info(self._ClassName(), info_msg)
        if self.echo_level > 2:
            info_msg  = colors.bold("abs_norm") + " = " + str(abs_norm) + " | "
            info_msg += colors.bold("abs_tol")  + " = " + str(self.abs_tolerance) + " || "
            info_msg += colors.bold("rel_norm") + " = " + str(rel_norm) + " | "
            info_msg += colors.bold("rel_tol")  + " = " + str(self.rel_tolerance)
            cs_tools.cs_print_info(self._ClassName(), info_msg)

        return is_converged

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "abs_tolerance" : 1e-5,
            "rel_tolerance" : 1e-5,
            "label"         : ""
        }""")
        this_defaults.AddMissingParameters(super(RelativeNormPreviousResidualConvergenceCriteria, cls)._GetDefaultSettings())
        return this_defaults

