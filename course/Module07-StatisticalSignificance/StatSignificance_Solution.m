% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/temp/JIDT/JIDT-copy/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/temp/JIDT/JIDT-copy/demos/octave');

% 0. Load/prepare the data:
N = 1000; % Number of samples to use
S = 1000; % Number of surrogates to generate
source = randn(N, 1); % assign random normal data to source
% destination = randn(N, 1); % assign random normal data to destination
coupling = 0.05;
destination = coupling .* source + (1 - coupling) .* randn(N, 1); % couple the destination to the source

% 1. Construct the calculator:
% calc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1');
% 2. Set any properties to non-default values:
% No properties were set to non-default values
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.setObservations(source, destination);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();
% 6. Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
measDist = calc.computeSignificance(S);

fprintf('MI_%s(col_0 -> col_1) = %.4f nats (null: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)\n', ...
	extractBetween(calcName, 'continuous.', '.'), result, measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), measDist.pValue, S);
histogram(javaMatrixToOctave(measDist.distribution), 50);
[n, edges] = histcounts(javaMatrixToOctave(measDist.distribution), 50);
line([result, result], [0, max(n)], 'Color', '[0 1 0]', 'LineWidth', 2);  % Mark in our measured MI

% Now add a nice title to the plot
calcName = string(calc.getClass().getName());
title(sprintf('Surrogate distribution for %d samples,\ncoupling=%.2f, %d surrogates, estimator=%s', N, coupling, S, extractBetween(calcName, 'continuous.', '.')));
xlabel('MI');
ylabel('count(MI)');
