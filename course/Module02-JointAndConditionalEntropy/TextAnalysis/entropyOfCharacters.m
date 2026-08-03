str = fileread('Seinfeld-scripts-textOnly.txt');
processedStr = regexprep(str, '[!"#$%&''()\*,-./;<=>?\[\\\]_`{}~]', ''); % Remove punctuation characters
processedStr = regexprep(processedStr, '[0-9]', ''); % Removed digits
processedStr = replace(processedStr, newline, ' '); % Replace newline characters with spaces
processedStr = lower(processedStr); % Convert all upper case into lower case
unique(processedStr)

% Check that this is adding a path to your scripts correctly:
addpath('../../Module1-IntroToInfoTheory/MatlabCode/completed');

% Compute the entropy of individual characters:
[result, symbols, probabilities] = entropyempirical(processedStr);
fprintf('Entropy of individual characters: %.4f bits\n', result);

% Now compute the info content of individual characters:
characterInfoContents = infocontent(probabilities);

% To just dump them to the screen:
%  characterInfoContents
% To display more nicely:
for ix = 1:length(symbols)
    fprintf('Info content of %s is %.4f bits\n', symbols(ix), characterInfoContents(ix));
end

% Now compute joint entropies for two characters:
%  Need a matrix with first column being first character, and second column
%  being the second
characterPairSamples = [processedStr(1:end-1)',processedStr(2:end)'];
pairEntropy = jointentropyempirical(characterPairSamples);
fprintf('Entropy of characters pairs: %.4f bits\n', pairEntropy);
% Finally compute the conditional entropy of the second character given the
% first:
conditionalEntropy = conditionalentropyempirical(characterPairSamples(:,2), characterPairSamples(:,1));
fprintf('Conditional entropy of character given previous: %.4f bits\n', conditionalEntropy);
