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

# = Example 7 - Ensemble method with transfer entropy on continuous data using Kraskov estimators =

#  Calculation of transfer entropy (TE) by supplying an ensemble of samples from multiple time series.
# We use continuous-valued data using the Kraskov-estimator TE calculator here.

from jpype import *
import random
import math
import os

# Change location of jar to match yours (we assume script is called from demos/python):
jarLocation = os.path.join(os.getcwd(), "..", "..", "infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Generate some random normalised data.
numObservations = 1000
covariance=0.4
numTrials=10
kHistoryLength=1

# Create a TE calculator and run it:
teCalcClass = JPackage("infodynamics.measures.continuous.kraskov").TransferEntropyCalculatorKraskov
teCalc = teCalcClass()
teCalc.setProperty("k", "4") # Use Kraskov parameter K=4 for 4 nearest points
teCalc.initialise(kHistoryLength) # Use target history length of kHistoryLength (Schreiber k)
teCalc.startAddObservations()

for trial in range(0,numTrials):
	# Create a new trial, with destArray correlated to
	#  previous value of sourceArray:
	sourceArray = [random.normalvariate(0,1) for r in range(numObservations)]
	destArray = [0] + [sum(pair) for pair in zip([covariance*y for y in sourceArray[0:numObservations-1]], \
		[(1-covariance)*y for y in [random.normalvariate(0,1) for r in range(numObservations-1)]] ) ]
	
	# Add observations for this trial:
	print("Adding samples from trial %d ..." % trial)
	teCalc.addObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray))

# We've finished adding trials:
print("Finished adding trials")
teCalc.finaliseAddObservations()

# Compute the result:
print("Computing TE ...")
result = teCalc.computeAverageLocalOfObservations()
# Note that the calculation is a random variable (because the generated
#  data is a set of random variables) - the result will be of the order
#  of what we expect, but not exactly equal to it; in fact, there will
#  be some variance around it (smaller than example 4 since we have more samples).
print("TE result %.4f nats; expected to be close to %.4f nats for these correlated Gaussians " % \
	(result, math.log(1.0/(1-math.pow(covariance,2)))))

# And here's how to pull the local TEs out corresponding to each input time
# series under the ensemble method (i.e. for multiple trials).
localTEs=teCalc.computeLocalOfPreviousObservations()
localValuesPerTrial = teCalc.getSeparateNumObservations()  # Need to convert to int for indices later
startIndex = 0
for localValuesInThisTrial in localValuesPerTrial:
	endIndex = startIndex + localValuesInThisTrial - 1
	print("Local TEs for trial %d go from array index %d to %d" % (trial, startIndex, endIndex))
	print("  corresponding to time points %d:%d (indexed from 0) of that trial" % (kHistoryLength, numObservations-1))
	# Access the local TEs for this trial as:
	localTEForThisTrial = localTEs[startIndex:endIndex]
	# Now update the startIndex before we go to the next trial
	startIndex = endIndex + 1
# And make a sanity check that we've looked at all of the local values here:
print("We've looked at %d local values in total, matching the number of samples we have (%d)" % (startIndex, teCalc.getNumObservations()))

