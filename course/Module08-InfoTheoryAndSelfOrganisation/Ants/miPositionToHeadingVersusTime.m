% This sample analysis looks at how the information in x position about
%  the heading of each ant changes across time.
% You can see that the MI increases when the ants are making trails,
% with the largest 3 increases corresponding to when the trails to each food
% source are being sustained.
% Subseuqent analysis investigates this in more subtle ways, but does not
% give us additional insights.

% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
positionsx = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ants/positionsx.txt');
positionsy = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ants/positionsy.txt');
headings = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ants/headings.txt');

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
% 2. Set any properties to non-default values:
% No properties were set to non-default values

misVersusTime = zeros(size(positionsx, 1), 1);

% Compute for samples taken across all pairs at each time step:
for t = 1:size(positionsx, 1)
    source = octaveToJavaDoubleArray(positionsx(t, :));
    destination = octaveToJavaDoubleArray(headings(t, :));

    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    % 4. Supply the sample data:
    calc.setObservations(source, destination);
    % 5. Compute the estimate:
    misVersusTime(t) = calc.computeAverageLocalOfObservations();

    fprintf('MI_Kraskov (KSG) alg. 2(time=%d) = %.4f nats\n', ...
        t, misVersusTime(t));
end

plot(misVersusTime);
xlabel('time');
ylabel('MI from x position to heading (nats)');
% Compare the times of the peaks in MI to when the trails are formed and collapse.
% What is this detecting?

% Some additional analysis to address more subtle questions follows:

% What if we put all of the positions versus headings into a single calculation, then pulled out averages for each time step?
% Also remove the first 125 time steps, since not all ants are out yet
positionsx = positionsx(126:end, :);
positionsy = positionsy(126:end, :);
headings = headings(126:end, :);
source = octaveToJavaDoubleArray(reshape(positionsx, size(positionsx, 1) * size(positionsx, 2), 1));
destination = octaveToJavaDoubleArray(reshape(headings, size(headings, 1) * size(headings, 2), 1));

calc.initialise();
calc.setObservations(source, destination);
result = calc.computeAverageLocalOfObservations();
localValues = calc.computeLocalOfPreviousObservations();
averageAtEachTime = mean(reshape(localValues, size(positionsx, 1), size(positionsx, 2)), 2); % Take row means for each time
figure();
plot(averageAtEachTime);
xlabel('time')
ylabel('Local MI (nats)')
title('MI for x positon to heading, averaged across all ants at each time step');
% This doesn't show much - looks like x position alone doesn't tell us a
% whole lot, when we put all data into a single calculation

% What if we bring x and y position together?
source = octaveToJavaDoubleMatrix([reshape(positionsx, size(positionsx, 1) * size(positionsx, 2), 1), reshape(positionsy, size(positionsy, 1) * size(positionsy, 2), 1)]);
destination = octaveToJavaDoubleMatrix(reshape(headings, size(headings, 1) * size(headings, 2), 1));

calc.initialise(2, 1);
calc.setObservations(source, destination);
result = calc.computeAverageLocalOfObservations();
localValues = calc.computeLocalOfPreviousObservations();
averageAtEachTime = mean(reshape(localValues, size(positionsx, 1), size(positionsx, 2)), 2); % Take row means for each time
figure();
plot(averageAtEachTime);
% This doesn't show much - looks like x position alone doesn't tell us a whole lot.


