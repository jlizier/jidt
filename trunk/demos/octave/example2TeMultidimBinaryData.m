% = Example 2 - Transfer entropy on multidimensional binary data =

% Simple transfer entropy (TE) calculation on multidimensional binary data using the discrete TE calculator.

% This example is important for Octave users, because it shows how to handle multidimensional arrays from Octave to Java (this is not as simple as single dimensional arrays in example 1 - it requires using supplied scripts to convert the array).

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Create many columns in a multidimensional array,
%  where the next time step (row 2) copies the value of the column on the left
%  from the previous time step (row 1):
twoDTimeSeriesOctave       = (rand(1, 100)>0.5)*1;
twoDTimeSeriesOctave(2, :) = [twoDTimeSeriesOctave(1,100), twoDTimeSeriesOctave(1, 1:99)];

% Things get a little tricky if we want to pass 2D arrays into Java.
% Unlike native Octave 1D arrays in Example 1, 
%  native Octave 2D+ arrays do not seem to get directly converted to java arrays,
%  so we use the supplied scripts to make the conversion (via org.octave.Matrix class in octave)
% Matlab handles the conversion automatically, so in Matlab this script just returns
%  the array that was passed in.
twoDTimeSeriesJavaInt = octaveToJavaIntMatrix(twoDTimeSeriesOctave);

% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.discrete.ApparentTransferEntropyCalculator', 2, 1);
teCalc.initialise();
% Add observations of transfer across one cell to the right per time step:
teCalc.addObservations(twoDTimeSeriesJavaInt, 1);
fprintf('The result should be close to 1 bit here, since we are executing copy operations of what is effectively a random bit to each cell here: ');
result2D = teCalc.computeAverageLocalOfObservations()

