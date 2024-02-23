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

b = 27; % 27 characters including space

for k = 1:10

    % adapting code from week 4 conditionalMIAsFunctionOfLag.m:
    historyCharacters = [];
    nextChar = data(k+1:end)'; % Make it a column vector
    for lag=1:k
        laggedCharSamples = data(k+1-lag:end-lag)'; % Make it a column vector
        historyCharacters = [historyCharacters, laggedCharSamples];
    end
    
    % Now we can symbolize the historyCharacters in the same way as jointentropy solution code:
    [histSymbolSet,~,combinedHistorySymbol] = unique(historyCharacters, 'rows');

    % And JIDT requires the symbols to start from 0, so adjust here:
    histSymbolSet = histSymbolSet - 1;
    combinedHistorySymbol = combinedHistorySymbol - 1;

    fprintf('Symbolisation of the length-%d history vectors results in %d unique symbols (versus %d state space size)\n', ...
        k, length(histSymbolSet), b^k);

    % Now we will computed MI between combinedHistorySymbol and nextChar

    % 1. Construct the calculator:
    calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', b, length(histSymbolSet), 0);
    % 2. No other properties to set for discrete calculators.
    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    % 4. Supply the sample data:
    calc.addObservations(nextChar, combinedHistorySymbol);
    % 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations();
    results(k) = result;
    % 6. Compute the (statistical significance via) null distribution analytically:
    % **NOTE** the analytic null distribution here will look different to
    % the basic AIS estimator that we are comparing to -- this is because
    % it's computed from the actual (smaller) state space size here (so is
    % more accurate)
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

