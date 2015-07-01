#!/bin/bash

# The demo is intended to be run from the parent directory
cd ..

# Make sure the latest example source file is compiled.
javac -classpath "../../infodynamics.jar" "infodynamics/demos/autoanalysis/AutoAnalyser.java"

# Run the example:
java -classpath ".:../../infodynamics.jar" infodynamics.demos.autoanalysis.AutoAnalyser

cd $OLDPWD

