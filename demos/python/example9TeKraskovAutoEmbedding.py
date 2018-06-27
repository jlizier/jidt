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

# = Example 9 - Transfer entropy on continuous data using Kraskov estimators with auto-embedding =

# Transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator,
# with automatic selection of embedding parameters 

from jpype import *
import random
import math
import numpy
import readFloatsFile
import os

# Change location of jar to match yours (we assume script is called from demos/python):
jarLocation = os.path.join(os.getcwd(), "..", "..", "infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Examine the heart-breath interaction that Schreiber originally looked at:
datafile = '../data/SFI-heartRate_breathVol_bloodOx.txt'
data = readFloatsFile.readFloatsFile(datafile)
# As numpy array:
A = numpy.array(data)
# Select data points 2350:3550, pulling out the relevant columns:
breathRate = A[2350:3551,1]; 
heartRate = A[2350:3551,0];

# Create a Kraskov TE calculator:
teCalcClass = JPackage("infodynamics.measures.continuous.kraskov").TransferEntropyCalculatorKraskov
teCalc = teCalcClass()

# Set properties for auto-embedding of both source and destination
#  using the Ragwitz criteria:
#  a. Auto-embedding method
teCalc.setProperty(teCalcClass.PROP_AUTO_EMBED_METHOD,
		teCalcClass.AUTO_EMBED_METHOD_RAGWITZ)
#  b. Search range for embedding dimension (k) and delay (tau)
teCalc.setProperty(teCalcClass.PROP_K_SEARCH_MAX, "6")
teCalc.setProperty(teCalcClass.PROP_TAU_SEARCH_MAX, "6")
# Since we're auto-embedding, no need to supply k, l, k_tau, l_tau here:
teCalc.initialise()
# Compute TE from breath (column 1) to heart (column 0) 
teCalc.setObservations(breathRate, heartRate)
teBreathToHeart = teCalc.computeAverageLocalOfObservations()

# Check the auto-selected parameters and print out the result:
optimisedK = int(teCalc.getProperty(teCalcClass.K_PROP_NAME))
optimisedKTau = int(teCalc.getProperty(teCalcClass.K_TAU_PROP_NAME))
optimisedL = int(teCalc.getProperty(teCalcClass.L_PROP_NAME))
optimisedLTau = int(teCalc.getProperty(teCalcClass.L_TAU_PROP_NAME))
print(("TE(breath->heart) was %.3f nats for (heart embedding:) k=%d," + \
		"k_tau=%d, (breath embedding:) l=%d,l_tau=%d optimised via Ragwitz criteria") % \
		(teBreathToHeart, optimisedK, optimisedKTau, optimisedL, optimisedLTau))

# Next, embed the destination only using the Ragwitz criteria:
teCalc.setProperty(teCalcClass.PROP_AUTO_EMBED_METHOD,
		teCalcClass.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)
teCalc.setProperty(teCalcClass.PROP_K_SEARCH_MAX, "6")
teCalc.setProperty(teCalcClass.PROP_TAU_SEARCH_MAX, "6")
# Since we're only auto-embedding the destination, we supply
#  source embedding here (to overwrite the auto embeddings from above):
teCalc.setProperty(teCalcClass.L_PROP_NAME, "1")
teCalc.setProperty(teCalcClass.L_TAU_PROP_NAME, "1")
# Since we're auto-embedding, no need to supply k and k_tau here:
teCalc.initialise()
# Compute TE from breath (column 1) to heart (column 0) 
teCalc.setObservations(breathRate, heartRate)
teBreathToHeartDestEmbedding = teCalc.computeAverageLocalOfObservations()

# Check the auto-selected parameters and print out the result:
optimisedK = int(teCalc.getProperty(teCalcClass.K_PROP_NAME))
optimisedKTau = int(teCalc.getProperty(teCalcClass.K_TAU_PROP_NAME))
print(("TE(breath->heart) was %.3f nats for (heart embedding:) k=%d," + \
		"k_tau=%d, optimised via Ragwitz criteria, plus (breath embedding:) l=1,l_tau=1") % \
		(teBreathToHeartDestEmbedding, optimisedK, optimisedKTau))

