# = Example 6 - Mutual information calculation with dynamic specification of calculator =

# This example shows how to write Python code to take advantage of the
# common interfaces defined for various information-theoretic calculators.
# Here, we use the common form of the infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  interface (which is never named here) to write common code into which we can plug
#  one of three concrete implementations (kernel estimator, Kraskov estimator or
#  linear-Gaussian estimator) by dynamically supplying the class name of
#  the concrete implementation.
#
# This is the Python equivalent to the demos/java/lateBindingDemo

from jpype import *
import random
import string
import numpy

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

#---------------------
# 1. Properties for the calculation (these are dynamically changeable):
# The name of the data file (relative to this directory)
datafile = '../data/4ColsPairedNoisyDependence-1.txt'
# List of column numbers for variables 1 and 2:
#  (you can select any columns you wish to be contained in each variable)
variable1Columns = [0,1] # array indices start from 0 in python
variable2Columns = [2,3]
# The name of the concrete implementation of the interface 
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  which we wish to use for the calculation.
# Note that one could use any of the following calculators (try them all!):
#  implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1" # MI([0,1], [2,3]) = 0.35507
#  implementingClass = "infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel"
#  implementingClass = "infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian"
implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1"

#---------------------
# 2. Load in the data (space separate numbers, one time step per line, each column is a variable)
f = open(datafile)
data = []
for line in f:
	data.append([float(x) for x in line.split()])
	
# As numpy array:
A = numpy.array(data)

# Pull out the columns from the data set which correspond to each of variable 1 and 2:
variable1 = A[:,variable1Columns]
variable2 = A[:,variable2Columns]

#--------------------
# 3. Dynamically instantiate an object of the given class:
# (in fact, all java object creation in python is dynamic - it has to be,
#  since the languages are interpreted. This makes our life slightly easier at this
#  point than it is in demos/java/lateBindingDemo where we have to handle this manually)
indexOfLastDot = string.rfind(implementingClass, ".")
implementingPackage = implementingClass[:indexOfLastDot]
implementingBaseName = implementingClass[indexOfLastDot+1:]
miCalcClass = eval('JPackage(\'%s\').%s' % (implementingPackage, implementingBaseName))
miCalc = miCalcClass()

#--------------------
# 4. Start using the MI calculator, paying attention to only
#  call common methods defined in the interface type
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  not methods only defined in a given implementation class.
# a. Initialise the calculator to use the required number of
#   dimensions for each variable:
miCalc.initialise(len(variable1Columns), len(variable2Columns))
# b. Supply the observations to compute the PDFs from:
miCalc.setObservations(variable1, variable2)
# c. Make the MI calculation:
miValue = miCalc.computeAverageLocalOfObservations()

print("MI calculator %s computed the joint MI as %.5f\n" % (implementingClass, miValue))

