% assumes processedStr holds the text as previously processed - you can run
%  the previous solution code entropyOfCharacters.m to pull this up

% Compute conditional MI as a function of lag:

maxLag = 5;

condMisVsLag = zeros(1, maxLag);
conditionalCharacters = [];

nextChar = processedStr(maxLag+1:end);
for lag=1:maxLag
    sourceChar = processedStr(maxLag+1-lag:end-lag);
    if (lag == 1)
        % Nothing to condition on, just compute MI
        condMisVsLag(lag) = mutualinformationempirical(sourceChar, nextChar);
    else
        condMisVsLag(lag) = conditionalmutualinformationempirical(sourceChar, nextChar, conditionalCharacters);
    end
    fprintf('cond MI over lag %d is %.4f\n', lag, condMisVsLag(lag));
    conditionalCharacters = [conditionalCharacters, sourceChar'];
end

figure();
plot(1:maxLag, condMisVsLag, 'rx');
xlabel('Lag')
ylabel('MI (bits)');
title('Average MI between characters separated by the given lag conditioned on intervening chars');

