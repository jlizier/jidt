% Copyright (C) 2023, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week10/ca54.txt');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 4, 1, 1, 1, 1);
% 2. No other properties to set for discrete calculators.

networkSize = 100;
results = zeros(networkSize,networkSize);
pValues = zeros(networkSize,networkSize);
% Compute for all pairs:
for d = 1:networkSize
    for s = 1:networkSize
		% For each source-dest pair:
		if (s == d)
            pValues(s,d) = 1;
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
        % 6. Compute the (statistical significance via) null distribution analytically:
        measDist = calc.computeSignificance();

        results(s,d) = result;
        pValues(s,d) = measDist.pValue;
	end
end

figure(); imagesc(results); xlabel('target'); ylabel('source'); h = colorbar; ylabel(h, 'Local TE (bits)'); title('TE between all cell pairs');
threshold = 0.1;
figure(); imagesc(results > threshold); xlabel('target'); ylabel('source'); h = colorbar; ylabel(h, sprintf('TE > %.2f', threshold)); title('TE between all cell pairs -- compared to threshold');
figure(); imagesc(data(1:100, 1:networkSize)); xlabel('target'); ylabel('time'); title('Raw CA values'); h = colorbar; ylabel(h, 'Cell value'); 
figure(); imagesc(pValues); xlabel('target'); ylabel('source'); title('p-value for TE between all cell pairs'); h = colorbar; ylabel(h, 'p-value');
caxis([0,0.05/(networkSize*(networkSize-1))]);