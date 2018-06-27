##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2012, Joseph T. Lizier
##  
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

# = Example 1 - Transfer entropy on binary data =

# Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator:

import jpype
import random
import numpy
import os

# Change location of jar to match yours (we assume script is called from demos/python):
jarLocation = os.path.join(os.getcwd(), "..", "..", "infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
jpype.startJVM(jpype.getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Generate some random binary data.
sourceArray = [random.randint(0,1) for r in range(100)]
destArray = [0] + sourceArray[0:99]
sourceArray2 = [random.randint(0,1) for r in range(100)]

# Create a TE calculator and run it:
teCalcClass = jpype.JPackage("infodynamics.measures.discrete").TransferEntropyCalculatorDiscrete
teCalc = teCalcClass(2,1)
teCalc.initialise()

# First use simple arrays of ints, which we can directly pass in:
teCalc.addObservations(sourceArray, destArray)
print("For copied source, result should be close to 1 bit : %.4f" % teCalc.computeAverageLocalOfObservations())
teCalc.initialise()
teCalc.addObservations(sourceArray2, destArray)
print("For random source, result should be close to 0 bits: %.4f" % teCalc.computeAverageLocalOfObservations())

# Next, demonstrate how to do this with a numpy array
teCalc.initialise()
# Create the numpy arrays:
sourceNumpy = numpy.array(sourceArray, dtype=numpy.int)
destNumpy = numpy.array(destArray, dtype=numpy.int)
# The above can be passed straight through to JIDT in python 2:
# teCalc.addObservations(sourceNumpy, destNumpy)
# But you need to do this in python 3:
sourceNumpyJArray = jpype.JArray(jpype.JInt, 1)(sourceNumpy.tolist())
destNumpyJArray = jpype.JArray(jpype.JInt, 1)(destNumpy.tolist())
teCalc.addObservations(sourceNumpyJArray, destNumpyJArray)
print("Using numpy array for copied source, result confirmed as: %.4f" % teCalc.computeAverageLocalOfObservations())

jpype.shutdownJVM() 

