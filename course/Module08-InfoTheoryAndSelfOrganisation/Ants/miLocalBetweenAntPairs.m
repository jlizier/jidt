% This sample analysis looks at how the information in x position
%  between all pairs of ants. We use pointwise values to 
%  rearrange this to look at how MI between the ants varies with time.
% You can see that the MI increases when the ants are making trails,
% with the largest 3 increases corresponding to when the trails to each food
% source are being sustained.

% Add JIDT jar library to the path
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ants/positionsx.txt');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 8, 8, 0);
% 2. No other properties to set for discrete calculators.
mUtils = javaObject('infodynamics.utils.MatrixUtils');

% Compute for all pairs:
avLocals = zeros(1, size(data, 1));
for s = 1:125
	for d = 1:125
		% For each source-dest pair:
		if (s == d)
			continue;
		end
		% Column indices start from 1 in Matlab:
		source = mUtils.discretise(octaveToJavaDoubleArray(data(:,s)), 8);
		destination = mUtils.discretise(octaveToJavaDoubleArray(data(:,d)), 8);

		% 3. Initialise the calculator for (re-)use:
		calc.initialise();
		% 4. Supply the sample data:
		calc.addObservations(source, destination);
		% 5. Compute the estimate:
		result = calc.computeAverageLocalOfObservations();
		% 6. Compute the (statistical significance via) null distribution analytically:
		% measDist = calc.computeSignificance();

		% fprintf('MI_Binned(col_%d -> col_%d) = %.4f bits (analytic p(surrogate > measured)=%.5f)\n', ...
		% 	s, d, result, measDist.pValue);
		locals = calc.computeLocalFromPreviousObservations(source, destination);
		avLocals = avLocals + javaMatrixToOctave(locals);
    end
    fprintf('Done source ant %d\n', s);
end
avLocals = avLocals ./ (125*124);
plot(avLocals);
xlabel('time');
ylabel('Average local MI between pairs (bits)');
