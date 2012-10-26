% function jIntArray = octaveToJavaIntArray(octaveArray)
%
% Convert a native octave array to a java int 1D array.
% 
% Assumes the JIDT jar is already on the java classpath - you will get a 
% java classpath error if this is not the case.
%
% Unfortunately octave-java doesn't seem to handle the conversion properly,
%  and we must convert the array from a double array first.

function jIntArray = octaveToJavaIntArray(octaveArray)

	% Convert to a java Double array first - it doesn't seem to work converting elements to integers directly
	jDoubleArray = octaveToJavaDoubleArray(octaveArray);

	% Then convert this double matrix to an integer matrix
	mUtils = javaObject('infodynamics.utils.MatrixUtils');
	jIntArray = mUtils.doubleToIntArray(jDoubleArray);

end

