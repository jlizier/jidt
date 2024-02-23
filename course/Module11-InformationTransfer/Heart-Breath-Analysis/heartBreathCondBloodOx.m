% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/data/SFI-heartRate_breathVol_bloodOx-extract.txt');
% Column indices start from 1 in Matlab:
source = octaveToJavaDoubleArray(data(:,1));
destination = octaveToJavaDoubleArray(data(:,2));
conditional = octaveToJavaDoubleArray(data(:,3));

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.ConditionalTransferEntropyCalculatorKraskov');
% 2. Set any properties to non-default values:
calc.setProperty('k_HISTORY', '2');
calc.setProperty('l_HISTORY', '2');
calc.setProperty('DELAY', '6');
calc.setProperty('COND_EMBED_LENGTHS', '2'); % Using the embedding length for optimal pairwise info transfer.
calc.setProperty('COND_TAUS', '3'); % Not relevant if we're only using default embedding of 1 on conditional
calc.setProperty('COND_DELAYS', '1');
calc.setProperty('DYN_CORR_EXCL', '21');
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.setObservations(source, destination, conditional);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();
% 6. Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
measDist = calc.computeSignificance(100);

fprintf('CTE_Kraskov (KSG)(heart -> breath | bloodOx) = %.4f nats (null: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)\n', ...
	result, measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), measDist.pValue, 100);

localCondTEs = calc.computeLocalOfPreviousObservations();

% First let's colour the transitions in breath rate with the local conditional TE values
delay = str2num(calc.getProperty('DELAY')); % Instead of hard-coding this in case of changes, we'll work it out
figure(); scatter(destination(delay:end-1), destination(1+delay:end), 10, localCondTEs(1+delay:end));
xlabel('breath\_rate(n)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local CTE (nats)');
title('Breath samples coloured by local cond TE from heart');
 
% Next let's make a scatter plot showing the heart rate values:
figure(); scatter(source(1:end-delay), destination(1+delay:end), 10, localCondTEs(1+delay:end));
xlabel('heart\_rate(n-lag+1)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local TE (nats)');
title(sprintf('Heart-breath samples (delay %d) coloured by local cond TE', delay));
