% Example 6 - Mutual information calculation with dynamic specification of calculator

% This example shows how to write Matlab/Octave code to take advantage of the
%  common interfaces defined for various information-theoretic calculators.
% Here, we use the common form of the infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  interface (which is never named here) to write common code into which we can plug
%  one of three concrete implementations (kernel estimator, Kraskov estimator or
%  linear-Gaussian estimator) by dynamically supplying the class name of
%  the concrete implementation.
%
% This is the Octave/Matlab equivalent to the demos/java/lateBindingDemo

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

%---------------------
% 1. Properties for the calculation (these are dynamically changeable):
% The name of the data file (relative to this directory)
datafile = '../data/4ColsPairedNoisyDependence-1.txt';
% List of column numbers for variables 1 and 2:
%  (you can select any columns you wish to be contained in each variable)
variable1Columns = [1,2]; % array indices start from 1 in octave/matlab
variable2Columns = [3,4];
% The name of the concrete implementation of the interface 
%  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  which we wish to use for the calculation.
% Note that one could use any of the following calculators (try them all!):
%  implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1";
%  implementingClass = "infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel";
%  implementingClass = "infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian";
implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1";

%---------------------
% 2. Load in the data
data = load(datafile);
% Pull out the columns from the data set which correspond to each of variable 1 and 2:
variable1 = data(:, variable1Columns);
variable2 = data(:, variable2Columns);

% 3. Dynamically instantiate an object of the given class:
% (in fact, all java object creation in octave/matlab is dynamic - it has to be,
%  since the languages are interpreted. This makes our life slightly easier at this
%  point than it is in demos/java/lateBindingDemo where we have to handle this manually)
miCalc = javaObject(implementingClass);

% 4. Start using the MI calculator, paying attention to only
%  call common methods defined in the interface type
%  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  not methods only defined in a given implementation class.
% a. Initialise the calculator to use the required number of
%   dimensions for each variable:
miCalc.initialise(length(variable1Columns), length(variable2Columns));
% b. Supply the observations to compute the PDFs from:
miCalc.setObservations(octaveToJavaDoubleMatrix(variable1), octaveToJavaDoubleMatrix(variable2));
% c. Make the MI calculation:
miValue = miCalc.computeAverageLocalOfObservations();

fprintf("MI calculator %s computed the joint MI as %.5f\n",
		implementingClass, miValue);

