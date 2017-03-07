from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import main_script, time

class Solution(main_script.Solution):
    
    def __init__(self):
        super().__init__()
        
    def Initialize(self):
        super().Initialize()
        
        self.InitializeBenchmarks()
        
    def InitializeBenchmarks(self):
        
        start = time.time()
        print("\nSeconds spent computing PerformanceAccessingGetYoung:")
        for i in range(0,10000):
            PerformanceUtilities().PerformanceAccessingGetYoung(self.spheres_model_part)        
        print(time.time() - start)
        
        start = time.time()
        print("Seconds spent computing PerformanceAccessingGetPropertiesYoungModulus:")
        for i in range(0,10000):
            PerformanceUtilities().PerformanceAccessingGetPropertiesYoungModulus(self.spheres_model_part)
        print(time.time() - start)
        
        start = time.time()
        print("Seconds spent computing PerformanceAccessingProcessInfo:")
        for i in range(0,10000):
            PerformanceUtilities().PerformanceAccessingProcessInfo(self.spheres_model_part)
        print(time.time() - start)
        print("\n")
        
Solution().Initialize()