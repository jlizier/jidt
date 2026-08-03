% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% Hard code our time series examples:
% variable = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
% variable = [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0];
variable = [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];
% variable = [0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1];

% Plot the data:
figure(1);
plot(variable, 'x', 'MarkerSize', 10); ylabel('x(n)', 'FontSize', 20); xlabel('n', 'FontSize', 20); axis([0,20,0,1.1]); grid;

% Set history length k
k = 3;

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete', 2, k);
% 2. No other properties to set for discrete calculators.
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.addObservations(variable);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();
% 6. Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
measDist = calc.computeSignificance(100);

fprintf('AIS_Discrete(col_0) = %.4f bits (null: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)\n', ...
	result, measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), measDist.pValue, 100);

% Pull out the local AIS values for each point in the time series:
localAISValues = calc.computeLocalFromPreviousObservations(variable);
figure(2);
% We only plot the local values from time index k onwards -- the AIS is undefined before this (localAISValues just fills these values with zeros)
plot(k+1:length(localAISValues), localAISValues(k+1:end), 'x', 'MarkerSize', 10); ylabel('AIS(n)', 'FontSize', 20); xlabel('n', 'FontSize', 20); axis([0,20,-0.5,2]); grid;

