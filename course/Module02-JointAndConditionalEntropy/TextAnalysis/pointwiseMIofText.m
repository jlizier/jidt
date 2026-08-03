% assumes processedStr holds the text as previously processed - you can run
%  the previous solution code entropyOfCharacters.m to pull this up

[result, jointSymbols, jointProbabilities, xSymbols, xProbabilities, ySymbols, yProbabilities] = mutualinformationempirical(processedStr(1:end-1), processedStr(2:end));

pointwiseMIs = zeros(length(xSymbols), length(ySymbols)); % Create array to store the pointwise MI values for each possible character pair
for firstCharIndex = 1:length(xSymbols)
    firstChar = xSymbols(firstCharIndex);
    probFirst = xProbabilities(firstCharIndex);
    for secondCharIndex = 1:length(ySymbols)
        secondChar = ySymbols(secondCharIndex);
        probSecond = yProbabilities(secondCharIndex);
        jointSymbolIndex = find((jointSymbols(:,1) == firstChar) & (jointSymbols(:,2) == secondChar));
        if isempty(jointSymbolIndex)
            pointwiseMIs(firstCharIndex, secondCharIndex) = 0; % No occurence, so set to 0
            continue;
        end
        probJoint = jointProbabilities(jointSymbolIndex);
        % Compute the pointwise MI from probJoint, probFirst and probSecond
        pointwiseMIs(firstCharIndex, secondCharIndex) = log2( probJoint ./ (probFirst .* probSecond) );
    end
end

figure();
imagesc(pointwiseMIs)
ylabel('First letter');
xlabel('Second letter');
h = colorbar;
h.Label.String = 'MI (bits)';
h.Label.Rotation = 90;
xticks(1:27)
xticklabels(ySymbols); % second letters - y - goes on x axis
h = gca();
h.XTickLabelRotation = 0; % Align the x labels properly
yticks(1:27)
yticklabels(xSymbols); % first letters - x - goes on y axis
title('MI between successive letters of text');
