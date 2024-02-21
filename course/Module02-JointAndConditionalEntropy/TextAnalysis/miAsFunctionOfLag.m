% assumes processedStr holds the text as previously processed - you can run
%  the previous solution code entropyOfCharacters.m to pull this up

% Compute MI as a function of lag:

maxLag = 10;

misVsLag = zeros(1, maxLag);

for lag=1:maxLag
    misVsLag(lag) = mutualinformationempirical(processedStr(1:end-lag), processedStr(1+lag:end));
end

plot(1:maxLag, misVsLag, 'rx');
xlabel('Lag')
ylabel('MI (bits)');
title('Average MI between characters separated by the given lag');
