% = Example 4 - Transfer entropy on continuous data using Kraskov estimators =

% Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
sourceArray=normrnd(0, 1, numObservations, 1);
destArray = [0; covariance*sourceArray(1:numObservations-1) + (1-covariance)*normrnd(0, 1, numObservations - 1, 1)];
sourceArray2=normrnd(0, 1, numObservations, 1); % Uncorrelated source
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
teCalc.initialise(1); % Use history length 1 (Schreiber k=1)
teCalc.setProperty('k', '4'); % Use Kraskov parameter K=4 for 4 nearest points
% Perform calculation with correlated source:
teCalc.setObservations(sourceArray, destArray);
result = teCalc.computeAverageLocalOfObservations();
% Note that the calculation is a random variable (because the generated
%  data is a set of random variables) - the result will be of the order
%  of what we expect, but not exactly equal to it; in fact, there will
%  be a large variance around it.
fprintf('TE result %.4f nats; expected to be close to %.4f nats for these correlated Gaussians\n', ...
    result, log(1/(1-covariance^2)));
% Perform calculation with uncorrelated source:
teCalc.initialise(); % Initialise leaving the parameters the same
teCalc.setObservations(sourceArray2, destArray);
result2 = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians\n', result2);


