% = Example 3 - Transfer entropy on continuous data using kernel estimators =

% Simple transfer entropy (TE) calculation on continuous-valued data using the (box) kernel-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
sourceArray=randn(numObservations, 1);
destArray = [0; covariance*sourceArray(1:numObservations-1) + (1-covariance)*randn(numObservations - 1, 1)];
sourceArray2=randn(numObservations, 1); % Uncorrelated source
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
teCalc.setProperty('NORMALISE', 'true'); % Normalise the individual variables
teCalc.initialise(1, 0.5); % Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
teCalc.setObservations(sourceArray, destArray);
% For copied source, should give something close to expected value for correlated Gaussians:
result = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f bits; expected to be close to %.4f bits for these correlated Gaussians but biased upwards\n', ...
    result, log(1/(1-covariance^2))/log(2));
teCalc.initialise(); % Initialise leaving the parameters the same
teCalc.setObservations(sourceArray2, destArray);
% For random source, it should give something close to 0 bits
result2 = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f bits; expected to be close to 0 bits for uncorrelated Gaussians but will be biased upwards\n', ...
    result2);

