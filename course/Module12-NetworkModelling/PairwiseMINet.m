% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week10/ca54.txt');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2, 2, 0);
% 2. No other properties to set for discrete calculators.

networkSize = 100;
results = zeros(networkSize,networkSize);
% Compute for all pairs:
for d = 1:networkSize
    for s = 1:networkSize
		% For each source-dest pair:
		if (s == d)
			continue;
		end
		% Column indices start from 1 in Matlab:
		source = octaveToJavaIntArray(data(:, s));
		destination = octaveToJavaIntArray(data(:, d));

		% 3. Initialise the calculator for (re-)use:
		calc.initialise();
		% 4. Supply the sample data:
		calc.addObservations(source, destination);
		% 5. Compute the estimate:
		result = calc.computeAverageLocalOfObservations();

        results(s,d) = result;
	end
end

figure(); imagesc(results); colorbar; xlabel('target'); ylabel('source'); title('MI between all cell pairs');
figure(); imagesc(data(1:100, 1:networkSize)); xlabel('target'); ylabel('time'); title('Raw CA values'); colorbar;
