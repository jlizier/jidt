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
calcClass = JPackage("infodynamics.measures.continuous.kraskov").MutualInfoCalculatorMultiVariateKraskov2
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
    results.append(result)

    print("MI_Kraskov (KSG) alg. 2(col_0 -> col_1, timeDiff=%d) = %.4f nats" %\
        (timeDiff, result))

# Now plot the results:
import matplotlib.pyplot as plt
plt.figure();
plt.scatter(timeDiffs, results, c='red', marker='x');
plt.title('MI (KSG) versus heart-to-breath time delay')
plt.xlabel('time delay')
plt.ylabel('MI (nats)')
# If we're running this in a script you will need the following
plt.show(block = False)

###############
# Now let's check the local values to see which samples contributed to the high MI.
###############
# First let's go back to zero lag:
calc.setProperty('TIME_DIFF', str(0));
calc.initialise();
calc.setObservations(source, destination);
localMIs = calc.computeLocalOfPreviousObservations(); # Computes the array of local values for each sample
# Now make a scatter plot of the data and the local MIs:
plt.figure(); plt.scatter(data[:,0], data[:,1], c=localMIs, marker='o', s=8);
plt.title('Heart-breath samples (lag 0) coloured by local MI')
plt.xlabel('Heart rate'); plt.ylabel('Breath rate'); plt.colorbar(label='Local MI (nats)')
# If we're running this in a script you will need the following
plt.show(block = False)

# Next check the local values for the delay which maximised MI:
maxIndex = numpy.argmax(numpy.array(results))
timeDiffForMax = timeDiffs[maxIndex];
calc.setProperty('TIME_DIFF', str(timeDiffForMax))
calc.initialise();
calc.setObservations(source, destination);
localMIs = calc.computeLocalOfPreviousObservations(); # Computes the array of local values for each sample
# Now make a scatter plot of the data and the local MIs:
plt.figure(); plt.scatter(data[:-timeDiffForMax,0], data[timeDiffForMax:,1], c=localMIs, marker='o', s=8);
plt.title('Heart-breath samples (lag %d) coloured by local MI' % timeDiffForMax)
plt.xlabel('Heart rate'); plt.ylabel('Breath rate'); plt.colorbar(label='Local MI (nats)')
# If we're running this in a script you will need the following
plt.show(block = False)

###############
# Finally, a demonstration of how a *multivariate* mutual information can be
# used here:
###############
lag1 = 0;
if timeDiffForMax != 0:
    lag2 = timeDiffForMax;
else:
    lag2 = 1;
mvCalc = calcClass()
# Initialise for calculation from 2 source variables to 1 target variable.
# In future this will be done by setting properties.
mvCalc.initialise(2,1);
# Set the observations, aligning to the maximum lag manually here since the lag is different for
# the two sources (lag1 = 0, lag2 is larger). Need to set the source as a matrix now.
mvSource = JArray(JDouble, 2)(numpy.column_stack((data[lag2:,0], data[:-lag2,0])).tolist()); # heart rate
laggedDestination = JArray(JDouble, 1)(data[lag2:,1].tolist()); # breath rate
mvCalc.setObservations(mvSource, laggedDestination);
result = mvCalc.computeAverageLocalOfObservations();
print('MI_KSG(heart(lags %d,%d) -> breath) = %.4f nats' %\
    (lag1, lag2, result));
# And let's take a look on a 3D scatter of how the points relate: (note the differences in how to make 3D scatter plots here!)
localMIs = mvCalc.computeLocalOfPreviousObservations();
fig = plt.figure();
ax = fig.add_subplot(projection='3d')
p = ax.scatter(data[lag2:,0], data[:-lag2,0], data[lag2:,1], c=localMIs, marker='o', s=8);
plt.title('Heart-breath samples (lags 0,%d) coloured by local multivariate MI' % timeDiffForMax)
ax.set_xlabel('Heart rate lag 0'); ax.set_ylabel('Heart rate lag %d' % lag2); ax.set_zlabel('Breath rate');
fig.colorbar(p, label='Local MI (nats)')
# If we're running this in a script you will need the following
plt.show(block = False)

# Now show all the plots and block on that:
plt.show()

