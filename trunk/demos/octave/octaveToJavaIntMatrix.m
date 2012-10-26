% function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)
%
% Convert a native octave matrix to a java int 2D array.
% 
% Assumes the JIDT jar is already on the java classpath - you will get a 
% java classpath error if this is not the case.
%
% Unfortunately octave-java doesn't seem to handle the conversion properly,
%  and we must convert the matrix from a double matrix first.

function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)

	% Convert to a java Double matrix first - it doesn't seem to work converting elements to integers directly
	jDoubleMatrix = octaveToJavaDoubleMatrix(octaveMatrix);

	% Then convert this double matrix to an integer matrix
	mUtils = javaObject('infodynamics.utils.MatrixUtils');
	jIntMatrix = mUtils.doubleToIntArray(jDoubleMatrix);

end

