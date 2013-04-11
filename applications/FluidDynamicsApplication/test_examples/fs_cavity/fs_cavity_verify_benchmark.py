import sys

kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)

import benchmarking

def Run():
  print "Running Fractional step element test: 2D cavity flow"
  return benchmarking.RunBenchmark("fs_benchmark.py", "fs_benchmark_ref.txt")

if __name__=='__main__':
  if Run():
    print "Test successful"
  else:
    print "Test failed"
