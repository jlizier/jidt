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

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
% 2. Set any properties to non-default values:
calc.setProperty('k_HISTORY', '2');
calc.setProperty('l_HISTORY', '2');

results = zeros(1,20);

for delay = 1:20
    calc.setProperty('DELAY', string(delay));
    calc.setProperty('DYN_CORR_EXCL', '15'); % Could set this before loop, but is fine here
    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    % 4. Supply the sample data:
    calc.setObservations(source, destination);
    % 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations();
    results(delay) = result;
    
    % 6. Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
    measDist = calc.computeSignificance(100);

    fprintf('TE_Kraskov (KSG)(col_0 -> col_1, delay %d) = %.4f nats (null: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)\n', ...
        delay, result, measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), measDist.pValue, 100);
end

figure(); plot(1:20, results, 'rx'); xlabel('delay (time units)'); ylabel('TE (nats)'); title('Heart to breath rate TE versus delay');

%delay = 6; % optimal value we determined.
[maxTE, delay] = max(results); % Find the maximum automatically
calc.setProperty('DELAY', string(delay));
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.setObservations(source, destination);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();

% Compute the local TE values
localTEs = calc.computeLocalOfPreviousObservations();

% First let's colour the transitions in breath rate with the TE values
figure(); scatter(destination(delay:end-1), destination(1+delay:end), 10, localTEs(1+delay:end));
xlabel('breath\_rate(n)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local TE (nats)');
title('Breath samples coloured by local TE from heart rate');

% Next let's make a scatter plot showing the heart rate values:
figure(); scatter(source(1:end-delay), destination(1+delay:end), 10, localTEs(1+delay:end));
xlabel('heart\_rate(n-delay+1)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local TE (nats)');
title(sprintf('Heart-breath samples (delay %d) coloured by local TE', delay));
