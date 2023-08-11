from jpype import *

# Add JIDT jar library to the path
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Construct an AutoAnalyser
calcClass = JPackage("infodynamics.demos.autoanalysis").AutoAnalyserLauncher
calc = calcClass(False)

