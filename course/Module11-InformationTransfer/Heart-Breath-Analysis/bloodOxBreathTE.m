% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/data/SFI-heartRate_breathVol_bloodOx-extract.txt');
% Column indices start from 1 in Matlab:
source = octaveToJavaDoubleArray(data(:,3));
destination = octaveToJavaDoubleArray(data(:,2));

% Plot the blood oxygen concentration time series
figure(); plot(data(:,3)); xlabel('n'); ylabel('Blood ox. conc.'); title('Blood oxygen concentration versus sample index');
% For the autocorrelation we will use autocorr:
figure(); autocorr(data(:,3), 'NumLags', 25); title('ACF for blood ox. conc. time series');
% And finally we plot it with the heart and breath rates
normedData = normalize(data); % To make the time series easier to compare
figure(); plot(normedData(:,1), 'r-'); hold on; plot(normedData(:,2), 'b-'); plot(normedData(:,3), 'g-'); hold off;
xlabel('n'); ylabel('Normalised time series values');
legend({'Heart rate', 'Breath rate', 'Blood ox. conc.'}); title('Normalised series versus sample index');

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
% 2. Set any properties to non-default values:
calc.setProperty('k_HISTORY', '2');
calc.setProperty('l_HISTORY', '2');
calc.setProperty('l_TAU', '3');

results = zeros(1,20);

for delay = 1:20
    calc.setProperty('DELAY', string(delay));
    calc.setProperty('DYN_CORR_EXCL', '21'); % Could set this before loop, but is fine here
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

figure(); plot(1:20, results, 'rx'); xlabel('delay'); ylabel('TE (nats)'); title('Blood ox to breath rate TE versus delay');

%delay = 1; % optimal value we determined.
[maxTE, delay] = max(results); % Find the maximum automatically
calc.setProperty('DELAY', string(delay));
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.setObservations(source, destination);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();

localTEs = calc.computeLocalOfPreviousObservations();

% First let's colour the transitions in breath rate with the TE values
figure(); scatter(destination(delay:end-1), destination(1+delay:end), 10, localTEs(1+delay:end));
xlabel('breath\_rate(n)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local TE (nats)');
title('Breath samples coloured by local TE from blood ox.');
 
% Next let's make a scatter plot showing the heart rate values:
figure(); scatter(source(1:end-delay), destination(1+delay:end), 10, localTEs(1+delay:end));
xlabel('blood\_ox(n-delay+1)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local TE (nats)');
title(sprintf('BloodOx-breath samples (delay %d) coloured by local TE', delay));
