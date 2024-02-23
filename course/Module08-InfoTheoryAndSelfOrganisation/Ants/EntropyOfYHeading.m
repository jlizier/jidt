% This sample analysis looks at how the diversity (measured by entropy)
%  in the positions of each ant changes across time.
% You can see that the diversity drops when the ants are making trails,
% with the largest 3 drops corresponding to when the trails to each food
% source are being sustained.

% Add JIDT jar library to the path
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ants/positionsy.txt');

% Transpose so that we have a column for each time step.
% This is different to our usual approach of each ant being a variable in a
% column.
data = data';

% 1. Construct the calculator, using 8 discrete bins
calc = javaObject('infodynamics.measures.discrete.EntropyCalculatorDiscrete', 8);
% 2. No other properties to set for discrete calculators.
mUtils = javaObject('infodynamics.utils.MatrixUtils');

% Compute for all time steps:
results = zeros(1,size(data,2));
for v = 1:size(data,2)
	% For each time step:
	% Column indices start from 1 in Matlab:
	variable = mUtils.discretise(octaveToJavaDoubleArray(data(:,v)), 8);

	% 3. Initialise the calculator for (re-)use:
	calc.initialise();
	% 4. Supply the sample data:
	calc.addObservations(variable);
	% 5. Compute the estimate:
	results(v) = calc.computeAverageLocalOfObservations();

	% fprintf('H_Binned(col_%d) = %.4f bits\n', ...
	%	v, result);
end

plot(results);
xlabel('time');
ylabel('Entropy of y-position (bits)');
