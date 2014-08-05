# = Example 1 - Transfer entropy on binary data =

# Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator:

from jpype import *
import random

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Generate some random binary data.
sourceArray = [random.randint(0,1) for r in xrange(100)]
destArray = [0] + sourceArray[0:99];
sourceArray2 = [random.randint(0,1) for r in xrange(100)]

# Create a TE calculator and run it:
teCalcClass = JPackage("infodynamics.measures.discrete").TransferEntropyCalculator
teCalc = teCalcClass(2,1)
teCalc.initialise()
# Since we have simple arrays of ints, we can directly pass these in:
teCalc.addObservations(destArray, sourceArray)
print("For copied source, result should be close to 1 bit : %.4f" % teCalc.computeAverageLocalOfObservations())
teCalc.initialise()
teCalc.addObservations(destArray, sourceArray2)
print("For random source, result should be close to 0 bits: %.4f" % teCalc.computeAverageLocalOfObservations())

shutdownJVM() 

