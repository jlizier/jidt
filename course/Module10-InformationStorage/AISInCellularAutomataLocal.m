% Copyright (C) 2019, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave/CellularAutomata');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week10/ca54.txt');

k=15;
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete', 2, k);
% 2. No other properties to set for discrete calculators.

% 3. Initialise the calculator for (re-)use:
calc.initialise();

% Compute for all variables:
for v = 1:10000
	% For each variable:
	% Column indices start from 1 in Matlab:
	variable = octaveToJavaIntArray(data(:, v));

	% 4. Supply the sample data:
	calc.addObservations(variable);
end

% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();
% 6. Compute the (statistical significance via) null distribution analytically:
measDist = calc.computeSignificance();
bias = measDist.getMeanOfDistribution();
	
fprintf('AIS_Discrete(all cols, k=%d) = %.4f bits from %d samples (bias %.4f, bias corrected %.4f)\n', ...
		k, result, calc.getNumObservations(), bias, result - bias);


% Now let's compute the local AIS values
locals = [];
for v = 1:10000
	% For each variable:
	% Column indices start from 1 in Matlab:
	variable = octaveToJavaIntArray(data(:, v));

	% 4. Compute the local AIS for this cell:
	localAISForThisCell = javaMatrixToOctave(calc.computeLocalFromPreviousObservations(variable));
    locals(:,v) = localAISForThisCell;
end
% Default plotting code:
% figure(2); imagesc(locals(1:100,1:100)); xlabel('cells'); ylabel('timestep'); colorbar;
% Plotting code with nicer colours:
figure(2); plotLocalInfoValues(locals(1:100,1:100)); xlabel('cells'); ylabel('timestep'); h = colorbar; ylabel(h, 'Local AIS (bits)');

