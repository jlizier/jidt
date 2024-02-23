% In this sample analysis we compute the MI from each pixel to the digit
% class:

% Add JIDT jar library to the path
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/MNIST/trialAndTest-processedData.mat');
originalClasses = data.trialAndTestData(:,1);
classes = originalClasses - 1; % Need classes to start from 0 for JIDT
numClasses = max(originalClasses);

% Binarise the pixel data:
threshold = 15; % From viewing histograms of pixel values this seems reasonable.
pixels = data.trialAndTestData(:,2:end) > threshold;

numPixels = size(pixels, 2);
numTrials = size(pixels, 1);
imageDimension = sqrt(numPixels); % Will be 28
fprintf('Loaded MNIST data with %d samples, classes %d:%d, and %d pixels per sample\n', length(originalClasses), min(originalClasses), max(originalClasses), size(pixels,2));
% Plot an example figure to check we've got the alignment correct:
figure();
imagesc(reshape(pixels(3,:), imageDimension, imageDimension)'); % Plotting the third digit as a sample

mis = zeros(numPixels, 1);

for p = 1:numPixels
    % 1. Construct the calculator: MI from binary pixel to class (alphabet size 10)
    calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2, 10, 0);
    % 2. No other properties to set for discrete calculators.
    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    % 4. Supply the sample data:
    calc.addObservations(pixels(:,p), classes);
    % 5. Compute the estimate:
    mis(p) = calc.computeAverageLocalOfObservations();
end

% Plot the MI from each pixel to the digit class:
figure();
imagesc(reshape(mis, imageDimension, imageDimension)');
h = colorbar;
title('MI from each pixel to the digit class');
xlabel('Pixel x index');
ylabel('Pixel y index');
ylabel(h, 'MI (bits)');
