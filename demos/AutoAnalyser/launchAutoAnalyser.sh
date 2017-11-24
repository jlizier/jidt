#!/bin/bash

# Need to change directory upwards since the launcher expects to run from JIDT root
cd ../..

# Run the AutoAnalyser Launcher:
java -classpath "infodynamics.jar" infodynamics.demos.autoanalysis.AutoAnalyserLauncher

