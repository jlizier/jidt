%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2016 Joseph T. Lizier
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

% = Example 10 - Conditional Transfer entropy on continuous multivariate data using Kraskov estimators =

% Conditional Transfer entropy (TE) calculation on multivariate continuous-valued data using the Kraskov-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar')

% Generate some random normalised data.
numObservations = 100000;
% Keep the sum of squares of these covariances below 1 to allow proper calculation of noise term.
covarianceToSource = 0.4;
covarianceToConds = [0.3,0.3];
if (sumsq([covarianceToSource, covarianceToConds]) >= 1)
	error('Sum of squares of the covariances must be < 1 here');
end
noiseCovar = sqrt(1 - sumsq([covarianceToSource, covarianceToConds]));

% Generate the random variables
sourceArray = randn(numObservations, 1);
condArray = randn(numObservations, length(covarianceToConds));
destArray = [0; covarianceToSource*sourceArray(1:numObservations-1,:) + covarianceToConds(1)*condArray(1:numObservations-1,1) + covarianceToConds(2)*condArray(1:numObservations-1,2) + noiseCovar*randn(numObservations-1, 1)];

% Expected results:
expectedConditional = -0.5 * log(noiseCovar .* noiseCovar ./ (1 - sumsq(covarianceToConds)));
expectedPairwise = -0.5 * log(1-covarianceToSource.^2);

% Create a conditional TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kraskov.ConditionalTransferEntropyCalculatorKraskov');
teCalc.initialise(1,1, ... % Destination embedding length (Schreiber k=1) and delays
		1,1, ...   % Source embedding length (Schreiber l=1) and delays
		1, ...     % Source-destination delay of 1 (default)
		octaveToJavaDoubleArray([1,1]), ... % Embedding lengths for each conditional variable
		octaveToJavaDoubleArray([1,1]), ... % Embedding delays for each conditional variable
		octaveToJavaDoubleArray([1,1]) ... % conditional-destination delays for each conditional variable
);
teCalc.setObservations(octaveToJavaDoubleArray(sourceArray), ...
	octaveToJavaDoubleArray(destArray), ...
	octaveToJavaDoubleMatrix(condArray));
% Perform calculation with correlated source, but no conditioning on other sources:
conditionalResult = teCalc.computeAverageLocalOfObservations();

% Create a pairwise TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
teCalc.initialise(); % Use default embeddings of 1 (e.g. Schreiber k=1 and l=1)
teCalc.setObservations(octaveToJavaDoubleArray(sourceArray), octaveToJavaDoubleArray(destArray));
% Perform calculation with correlated source, but no conditioning on other sources:
pairwiseResult = teCalc.computeAverageLocalOfObservations();

% Note that the calculation is a random variable (because the generated
%  data is a set of random variables) - the result will be of the order
%  of what we expect, but not exactly equal to it; in fact, there will
%  be some variance around it. It will probably be biased down here
%  due to small correlations between the supposedly uncorrelated variables.
fprintf('From %d samples:\nTE result conditional result = %.4f nats, pairwise = %.4f nats;\nexpected around %.4f nats (conditional) and %.4f nats (pairwise) for the correlated Gaussians\n', ...
	numObservations, conditionalResult, pairwiseResult, expectedConditional, expectedPairwise);

clear teCalc

