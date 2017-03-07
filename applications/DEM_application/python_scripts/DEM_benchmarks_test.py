from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import main_script_test, time, os, sys

import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi

sys.path.insert(0,'')
# DEM Application

benchmark_number = int(sys.argv[1])

listDISCONT   = list(range(1,12))
listROLLFR    = list(range(12,13))
listDEMFEM    = list(range(13,18))
listCONT      = list(range(20,27))
listDISclZHAO = [30,32]
listDISclRK   = [31,33]

if benchmark_number in listDISCONT:
    import DEM_explicit_solver_var as DEM_parameters
elif benchmark_number in listROLLFR:
    import DEM_explicit_solver_var_ROLLFR as DEM_parameters
elif benchmark_number in listDEMFEM:
    import DEM_explicit_solver_var_DEMFEM as DEM_parameters
elif benchmark_number in listCONT:
    import DEM_explicit_solver_var_CONT as DEM_parameters
elif benchmark_number == 27:
    import DEM_explicit_solver_var_UCS as DEM_parameters
elif benchmark_number == 28:
    import DEM_explicit_solver_var_PENDULO3D as DEM_parameters
elif benchmark_number in listDISclZHAO:
    import DEM_explicit_solver_var_DISclZHAO as DEM_parameters
elif benchmark_number in listDISclRK:
    import DEM_explicit_solver_var_DISclRK as DEM_parameters
else:
    print('Benchmark number does not exist')
    sys.exit()

class Solution(main_script_test.Solution):
    
    def __init__(self):
        super().__init__()
        
    def Initialize(self):
        super().Initialize()
        
        self.SetInitialData(coeff_of_restitution_iteration)
        
    def SetInitialData(self, coeff_of_restitution_iteration):
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration)
        '''print(coeff_of_restitution_iteration)
        print("ese")'''
        
    def GetMpFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM"
    
    def GetInletFilename(self):
        return 'benchmark' + "DEM_Inlet"
    
    def BeforeSolveOperations(self, time):
        '''print(time)
        print(dt)
        print('check')'''
        benchmark.ApplyNodalRotation(self.spheres_model_part, time, dt)
        
    def BeforePrintingOperations(self):
        #print(dt)
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.output_time_step, dt)
        
    def BeforeFinalizeOperations(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)
    
    def InitializeTimeStep(self):
        return dt
    
    def PreviousOperations(self):
        DEM_parameters.FinalTime, dt, self.output_time_step, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
        '''print(number_of_points_in_the_graphic)
        print(number_of_coeffs_of_restitution)
        print(self.output_time_step)
        oiga'''
nodeplotter = 0  
number_of_coeffs_of_restitution = 0

if benchmark_number in listDISCONT:
    nodeplotter = 1  

solution = Solution()

solution.PreviousOperations()
DEM_parameters.FinalTime, dt, _, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
'''print(DEM_parameters.FinalTime)
print(dt)
print(number_of_points_in_the_graphic)
print(number_of_coeffs_of_restitution)'''

benchmark_class_name = 'Benchmark' + str(benchmark_number)
benchmark_class = getattr(DBC, benchmark_class_name)
benchmark = benchmark_class()

for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    
    for iteration in range(1, number_of_points_in_the_graphic + 1):
    
        solution.Initialize()
        solution.SetInitialData(coeff_of_restitution_iteration)
        solution.RunMainTemporalLoop()
        solution.Finalize()
        print("hola")
        
    print("adios")
    print(number_of_points_in_the_graphic)
    print(dt)
    #colti
    benchmark.print_results(number_of_points_in_the_graphic, dt)
