% = Example 5 - Multivariate transfer entropy on binary data =

% Multivariate transfer entropy (TE) calculation on binary data using the discrete TE calculator:

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random binary data.
% Note that we need the *1 to make this a number not a Boolean,
%  otherwise this will not work (as it cannot match the method signature)
numObservations = 100;
sourceArray=(rand(numObservations,2)>0.5)*1;
sourceArray2=(rand(numObservations,2)>0.5)*1;
% Destination variable takes a copy of the first bit of the source in bit 1,
%  and an XOR of the two bits of the source in bit 2:
destArray = [0, 0; sourceArray(1:numObservations-1, 1), xor(sourceArray(1:numObservations-1, 1), sourceArray(1:numObservations-1, 2))];
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculator', 4, 1);
teCalc.initialise();
% We need to construct the joint values of the dest and source before we pass them in,
%  and need to use the matrix conversion routine when calling from Matlab/Octave:
mUtils= javaObject('infodynamics.utils.MatrixUtils');
teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaIntMatrix(sourceArray), 2), ...
		mUtils.computeCombinedValues(octaveToJavaIntMatrix(destArray), 2));
fprintf('For source which the 2 bits are determined from, result should be close to 2 bits : ');
result = teCalc.computeAverageLocalOfObservations()
teCalc.initialise();
teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaIntMatrix(sourceArray2), 2), ...
		mUtils.computeCombinedValues(octaveToJavaIntMatrix(destArray), 2));
fprintf('For random source, result should be close to 0 bits in theory: ');
result2 = teCalc.computeAverageLocalOfObservations()
fprintf('\nThe result for random source is inflated towards 0.3 due to finite observation length (%d).\nOne can verify that the answer is consistent with that from a\nrandom source by checking: teCalc.computeSignificance(1000); ans.pValue\n', teCalc.getNumObservations());

