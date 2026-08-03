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
calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
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

	fprintf('MI_KSG(col_0 -> col_1) = %.4f nats\n', ...
		result);
end

figure();
plot(timeDiffs, results, 'x');
title("MI (KSG) versus heart-to-breath time delay");
xlabel("time delay");
ylabel("MI (nats)");

%%%%%%%%%%%%%%%
% Now let's check the local values to see which samples contributed to the high MI.
%%%%%%%%%%%%%%%
% First let's go back to zero lag:
calc.setProperty('TIME_DIFF', string(0));
calc.initialise();
calc.setObservations(source, destination);
localMIs = calc.computeLocalOfPreviousObservations();
figure();
scatter(data(:,1), data(:,2), 10, localMIs);
h = colorbar;
xlabel('Heart rate'); ylabel('Breath rate'); ylabel(h, 'Local MI (nats)');
title('Heart-breath samples (lag 0) coloured by local MI');

% Next check the local values for the delay which maximised MI:
[~,maxIndex] = max(results);
timeDiffForMax = timeDiffs(maxIndex);
calc.setProperty('TIME_DIFF', string(timeDiffForMax));
calc.initialise();
calc.setObservations(source, destination);
localMIs = calc.computeLocalOfPreviousObservations();
figure();
scatter(data(1:end-timeDiffForMax,1), data(1+timeDiffForMax:end,2), 10, localMIs);
h = colorbar;
xlabel('Heart rate'); ylabel('Breath rate'); ylabel(h, 'Local MI (nats)');
title(sprintf('Heart-breath samples (lag %d) coloured by local MI', timeDiffForMax));

%%%%%%%%%%%%%%%
% Finally, a demonstration of how a *multivariate* mutual information can be
% used here:
%%%%%%%%%%%%%%%
lag1 = 0;
if timeDiffForMax ~= 0
    lag2 = timeDiffForMax;
else
    lag2 = 1;
end
mvCalc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
% Initialise for calculation from 2 source variables to 1 target variable.
% In future this will be done by setting properties.
mvCalc.initialise(2,1);
% Set the observations, aligning to the maximum lag manually here since the lag is different for
% the two sources (lag1 = 0, lag2 is larger). Need to set the source as a matrix now.
mvSource = octaveToJavaDoubleMatrix([data(1+lag2:end,1), data(1:end-lag2,1)]); % heart rate
laggedDestination = octaveToJavaDoubleArray(data(1+lag2:end,2)); % breath rate
mvCalc.setObservations(mvSource, laggedDestination);
result = mvCalc.computeAverageLocalOfObservations();
fprintf('MI_KSG(heart(lags %d,%d) -> breath) = %.4f nats\n', ...
    lag1, lag2, result);
% And let's take a look on a 3D scatter of how the points relate:
localMIs = mvCalc.computeLocalOfPreviousObservations();
figure(); scatter3(data(1+lag2:end,1), data(1:end-lag2,1), data(1+lag2:end,2), 20, localMIs);
xlabel('heart lag 0');
ylabel(sprintf('heart lag %d', lag2));
zlabel('breath');
