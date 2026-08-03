% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/temp/jidt-master/jidt/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/temp/jidt-master/jidt/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/temp/jidt-master/jidt/demos/octave/CellularAutomata/ca54.txt');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 15, 1, 1, 1, 1);
% 2. No other properties to set for discrete calculators.

% 3. Initialise the calculator for (re-)use:
calc.initialise();

% Compute for all pairs, with source one cell to left of target:
for d = 1:10000
	s = mod (d - 2, 10000) + 1;
	% For each source-dest pair:
	if (s == d)
		continue;
	end
	% Column indices start from 1 in Matlab:
	source = octaveToJavaIntArray(data(:, s));
	destination = octaveToJavaIntArray(data(:, d));

	% 4. Supply the sample data:
	calc.addObservations(source, destination);
end

% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();

fprintf('TE_Discrete(one cell to right) = %.4f bits from %d samples\n', ...
	result, calc.getNumObservations());

