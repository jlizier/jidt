#!/bin/bash

# Make sure the latest example source file is compiled:
#  Note: this does not nulify our claim that you can compile the code
#  once and then only change the properties file - we only compile here
#  every time this script is run so that we can capture any changes the user
#  made to the demo source code. The demo as written should be compiled once
#  then one can make dynamic changes to the props file and simply
#  run the class file without recompiling.
javac -classpath "../../infodynamics.jar" "infodynamics/demos/Example10GPUBenchmark.java"

# Run the example, feeding in the file names as the command line argument
file1=WhiteNoise.txt
file2=Correlated2D.txt
java -classpath ".:../../infodynamics.jar" infodynamics.demos.Example10GPUBenchmark $file1 $file2

# Plot the results
python plotExample10BenchmarkResults.py $file1 $file2
