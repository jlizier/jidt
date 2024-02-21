% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/temp/JIDT/JIDT-copy/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/temp/JIDT/JIDT-copy/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/temp/JIDT/JIDT-copy/demos/data/2coupledRandomCols-1.txt');
% Column indices start from 1 in Matlab:
source = octaveToJavaDoubleArray(data(:,1));
destination = octaveToJavaDoubleArray(data(:,2));

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
results = [];
kNNs = 4:15;
for kNN = kNNs

	% 2. Set any properties to non-default values:
	calc.setProperty('TIME_DIFF', '1');
	calc.setProperty('k', string(kNN));
	% 3. Initialise the calculator for (re-)use:
	calc.initialise();
	% 4. Supply the sample data:
	calc.setObservations(source, destination);
	% 5. Compute the estimate:
	result = calc.computeAverageLocalOfObservations();
	results(end+1) = result;

	fprintf('MI_Kraskov (KSG) alg. 2(col_0 -> col_1, k=%d) = %.4f nats\n', ...
		kNN, result);
end

plot(kNNs, results, 'x-');
title("MI versus k nearest neighbours");
xlabel("kNNs");
ylabel("MI (nats)");
axis([4, 16, 0, max(results)*1.1]);

