% Copyright (C) 2023, Joseph T. Lizier
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
figure(1); imagesc(data(1:100,1:100)); xlabel('cells'); ylabel('timestep'); h = colorbar; ylabel(h, 'Cell values');
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

	localTEForThisCell = calc.computeLocalFromPreviousObservations(source, destination);
	locals(:,d) = localTEForThisCell;
end
figure(2); plotLocalInfoValues(locals(1:100,1:100)); xlabel('cells'); ylabel('timestep'); h = colorbar; ylabel(h, 'Local TE (bits)');

