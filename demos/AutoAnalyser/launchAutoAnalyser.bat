@ECHO OFF

REM Need to change directory upwards since the launcher expects to run from JIDT root
cd ..\..

REM Run the AutoAnalyser launcher:
java -classpath "..\..\infodynamics.jar" infodynamics.demos.autoanalysis.AutoAnalyserLauncher

