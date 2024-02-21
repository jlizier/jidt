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
dataRaw = readFloatsFile.readFloatsFile("/home/joseph/JIDT/infodynamics-dist-1.6/demos/data/SFI-heartRate_breathVol_bloodOx-extract.txt")
# As numpy array:
data = numpy.array(dataRaw)
source = JArray(JDouble, 1)(data[:,0].tolist())
destination = JArray(JDouble, 1)(data[:,1].tolist())

# 1. Construct the calculator:
calcClass = JPackage("infodynamics.measures.continuous.gaussian").MutualInfoCalculatorMultiVariateGaussian
calc = calcClass()

results = []
timeDiffs = range(0,16)

for timeDiff in timeDiffs:
    # 2. Set any properties to non-default values:
    calc.setProperty("TIME_DIFF", str(timeDiff))
    # 3. Initialise the calculator for (re-)use:
    calc.initialise()
    # 4. Supply the sample data:
    calc.setObservations(source, destination)
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    results.append(result);

    print("MI_Gaussian(col_0 -> col_1, timeDiff=%d) = %.4f nats" %\
        (timeDiff, result))

# Now plot the results:
import matplotlib.pyplot as plt
plt.figure();
plt.scatter(timeDiffs, results, c='red', marker='x');
plt.title('MI (Gaussian) versus heart-to-breath time delay')
plt.xlabel('time delay')
plt.ylabel('MI (nats)')
# If we're running this in a script you will need the following
plt.show()

