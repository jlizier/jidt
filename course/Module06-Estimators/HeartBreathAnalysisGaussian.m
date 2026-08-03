% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/temp/JIDT/JIDT-copy/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/temp/JIDT/JIDT-copy/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/temp/JIDT/JIDT-copy/demos/data/SFI-heartRate_breathVol_bloodOx-extract.txt');
% Column indices start from 1 in Matlab:
source = octaveToJavaDoubleArray(data(:,1)); % heart rate
destination = octaveToJavaDoubleArray(data(:,2)); % breath rate

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
results = [];
timeDiffs = 0:15;
for timeDiff = timeDiffs
	% 2. Set any properties to non-default values:
	calc.setProperty('TIME_DIFF', string(timeDiff));
	% 3. Initialise the calculator for (re-)use:
	calc.initialise();
	% 4. Supply the sample data:
	calc.setObservations(source, destination);
	% 5. Compute the estimate:
	result = calc.computeAverageLocalOfObservations();
	results(end+1) = result;

	fprintf('MI_Gaussian(col_0 -> col_1, timeDiff=%d) = %.4f nats\n', ...
		timeDiff, result);
end

plot(timeDiffs, results, 'x-');
title("MI (Gaussian) versus heart-to-breath time delay");
xlabel("time delay");
ylabel("MI (nats)");

