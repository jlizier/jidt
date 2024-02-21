from jpype import *
import numpy
import sys
# Our python data file readers are a bit of a hack, python users will do better on this:
sys.path.append("/home/joseph/JIDT/infodynamics-dist-1.6/demos/python")
import readFloatsFile

# Add JIDT jar library to the path
jarLocation = "/home/joseph/JIDT/infodynamics-dist-1.6/infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# 0. Load/prepare the data:
dataRaw = readFloatsFile.readFloatsFile("/home/joseph/JIDT/infodynamics-dist-1.6/demos/data/2coupledRandomCols-1.txt")
# As numpy array:
data = numpy.array(dataRaw)
source = JArray(JDouble, 1)(data[:,0].tolist())
destination = JArray(JDouble, 1)(data[:,1].tolist())

# 1. Construct the calculator:
calcClass = JPackage("infodynamics.measures.continuous.kraskov").MutualInfoCalculatorMultiVariateKraskov2
calc = calcClass()

results = []
ks = range(4,16)

for k in ks:
    # 2. Set any properties to non-default values:
    calc.setProperty("TIME_DIFF", "1")
    calc.setProperty("k", str(k))
    # 3. Initialise the calculator for (re-)use:
    calc.initialise()
    # 4. Supply the sample data:
    calc.setObservations(source, destination)
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    results.append(result);

    print("MI_Kraskov (KSG) alg. 2(col_0 -> col_1, k=%d) = %.4f nats" %\
        (k, result))

# Now plot the results:
import matplotlib.pyplot as plt
plt.figure();
plt.scatter(ks, results, c='red', marker='x');
plt.title('MI versus k nearest neighbours')
plt.xlabel('kNNs')
plt.ylabel('MI (bits)')
plt.axis([4, 16, 0, max(results)*1.1])
# If we're running this in a script you will need the following
plt.show()

