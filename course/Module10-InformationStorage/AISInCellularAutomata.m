% Copyright (C) 2019, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week10/ca54.txt');

for k = 1:20
	% 1. Construct the calculator:
	calc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete', 2, k);
	% 2. No other properties to set for discrete calculators.

	% 3. Initialise the calculator for (re-)use:
	calc.initialise();

	% Compute for all variables:
	for v = 1:10000
		% For each variable:
		% Column indices start from 1 in Matlab:
		variable = octaveToJavaIntArray(data(:, v));

		% 4. Supply the sample data:
		calc.addObservations(variable);
	end

	% 5. Compute the estimate:
	result = calc.computeAverageLocalOfObservations();
	results(k) = result;
	% 6. Compute the (statistical significance via) null distribution analytically:
	measDist = calc.computeSignificance();
	bias(k) = measDist.getMeanOfDistribution();
	
	fprintf('AIS_Discrete(all cols, k=%d) = %.4f bits from %d samples (bias %.4f, bias corrected %.4f)\n', ...
			k, result, calc.getNumObservations(), bias(k), result - bias(k));
end

biasCorrectedAIS = results - bias;
[maxAIS, optimalK] = max(biasCorrectedAIS);
fprintf('Optimal k=%d, giving bias corrected AIS(k=%d)=%.4f bits\n', optimalK,optimalK,biasCorrectedAIS(optimalK));

plot(results, 'rx');
hold on;
plot(bias, 'bo');
plot(biasCorrectedAIS, 'g+');
hold off;
xlabel('k');
ylabel('AIS (bits)');
legend('AIS raw', 'bias', 'bias corrected AIS');



