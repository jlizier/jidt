% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
str = fileread('../Module2-JointAndConditionalEntropy/TextAnalysis/Seinfeld-scripts-textOnly.txt'); % Load the data
processedStr = regexprep(str, '[!"#$%&''()\*,-./;<=>?\[\\\]_`{}~]', ''); % Remove punctuation characters
% processedStr = erasePunctuation(str); % Remove punctuation characters -- needs Text Analytics Toolbox
processedStr = regexprep(processedStr, '[0-9]', ''); % Removed digits
processedStr = replace(processedStr, newline, ' '); % Replace newline characters with spaces
processedStr = lower(processedStr); % Convert all upper case into lower case

data = processedStr - 'a'; % Converts the letters to an integer representing their offset from 'a'
% But spaces have been turned into the value ' ' - 'a' (i.e. the difference in character encodings of space and 'a').
% We need to replace this with the integer value 26, such that all of our symbols are consecutive for JIDT:
locationsOfSpace = find(data == ' ' - 'a');
data(locationsOfSpace) = 26;

for k = 1:4
    % 1. Construct the calculator:
    calc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete', 27, k);
    % 2. No other properties to set for discrete calculators.

    % 3. Initialise the calculator for (re-)use:
    calc.initialise();

    % 4. Supply the sample data:
    calc.addObservations(data);

    % 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations();
    results(k) = result;
    % 6. Compute the (statistical significance via) null distribution analytically:
    measDist = calc.computeSignificance();
    bias(k) = measDist.getMeanOfDistribution();
    
    fprintf('AIS_Discrete(k=%d) = %.4f bits from %d samples (bias %.4f, bias corrected %.4f)\n', ...
            k, result, calc.getNumObservations(), bias(k), result - bias(k));
end

plot(results, 'rx');
hold on;
plot(bias, 'bo');
biasCorrectedAIS = results - bias;
plot(biasCorrectedAIS, 'g+');
hold off;
xlabel('k');
ylabel('AIS (bits)');
legend('AIS raw', 'bias', 'bias corrected AIS');
title('AIS versus k for written English sample');


