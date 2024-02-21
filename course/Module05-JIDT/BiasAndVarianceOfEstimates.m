% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6/demos/octave');

numSamples = 10;
numCalcs = 50;
results = zeros(numCalcs, 1);

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2, 2, 0);

for r = 1 : numCalcs
    % 0. Load/prepare the data:
    % Column indices start from 1 in Matlab:
    source = randi(2,numSamples,1)-1;
    destination = source;

    % 2. No other properties to set for discrete calculators.
    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    % 4. Supply the sample data:
    calc.addObservations(source, destination);
    % 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations();
    results(r) = result;
    
    fprintf('MI_Discrete(col_0 -> col_1) = %.4f bits\n', ...
	    result);
end

fprintf('Results have mean %.3f bits (bias of %.3f bits) and variance %.3f bits\n', mean(results), mean(results) - 1, var(results));

