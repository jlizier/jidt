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
calc = javaObject('infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel');
results = [];
kernels = 0.1:0.05:1.0;
for kernelWidth = kernels

	% 2. Set any properties to non-default values:
	calc.setProperty('TIME_DIFF', '1');
	calc.setProperty('KERNEL_WIDTH', string(kernelWidth));
	% 3. Initialise the calculator for (re-)use:
	calc.initialise();
	% 4. Supply the sample data:
	calc.setObservations(source, destination);
	% 5. Compute the estimate:
	result = calc.computeAverageLocalOfObservations();
	results(end+1) = result;

	fprintf('MI_Kernel(col_0 -> col_1) = %.4f bits\n', ...
		result);
end

plot(kernels, results, 'x-');
title("MI versus kernel width");
xlabel("Kernel width");
ylabel("MI (bits)");

