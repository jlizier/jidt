% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/TeachingPlayground/CSYS5030/infodynamics-dist-1.5/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/TeachingPlayground/CSYS5030/infodynamics-dist-1.5/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/infodynamics-dist-1.5/demos/data/SFI-heartRate_breathVol_bloodOx-extract.txt');
% Column indices start from 1 in Matlab:
variable = octaveToJavaDoubleArray(data(:,2));

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov');
% 2. Set any properties to non-default values:
calc.setProperty('k_HISTORY', '2');
calc.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS');
calc.setProperty('AUTO_EMBED_K_SEARCH_MAX', '10');
calc.setProperty('AUTO_EMBED_TAU_SEARCH_MAX', '10');
calc.setProperty('DYN_CORR_EXCL', '15');
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.setObservations(variable);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();

fprintf('AIS_Kraskov (KSG)(col_1) = %.4f nats\n', ...
	result);

figure(); plot(1:length(variable), variable, 'r-'); xlabel('n'); ylabel('breath rate'); title('Breath rate versus sample index');

localAIS = calc.computeLocalOfPreviousObservations();
k_str = calc.getProperty('k_HISTORY');
k = str2num(k_str);
figure(); scatter(variable, localAIS); xlabel('breath rate'); ylabel('local AIS (nats)'); title('Local AIS versus breath rate');

figure(); scatter(1:length(variable), variable, 10, localAIS);
xlabel('n'); ylabel('breath rate'); h = colorbar; ylabel(h, 'Local AIS (nats)'); title('Breath rate versus sample index');
figure(); scatter(variable(1:end-1), variable(2:end), 10, localAIS(2:end));
xlabel('breath\_rate(n)'); ylabel('breath\_rate(n+1)'); h = colorbar; ylabel(h, 'Local AIS (nats)'); title('Breath rate next versus previous values');



