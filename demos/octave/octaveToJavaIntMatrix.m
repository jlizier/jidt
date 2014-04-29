% function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)
%
% Convert a native octave matrix to a java int 2D array.
% 
% Assumes the JIDT jar is already on the java classpath - you will get a 
% java classpath error if this is not the case.
%

function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		if ((rows(octaveMatrix)*columns(octaveMatrix)) > 1)
			% Do this the normal way
			tmp = javaObject('infodynamics.utils.OctaveMatrix');
			tmp.loadIntData(reshape(octaveMatrix,1,rows(octaveMatrix)*columns(octaveMatrix)),[rows(octaveMatrix), columns(octaveMatrix)]);
			jIntMatrix = tmp.asIntMatrix();
		else
			% For length 1 arrays, we need to perform a hack here or else
			%  java thinks the length-one array is a scalar.
			% See octaveToJavaDoubleArray for a further description.
			% So instead we'll do this the slow way (doesn't matter for one element only)
			jIntMatrix = javaArray('java.lang.Integer', 1, 1);
			jIntMatrix(1, 1) = octaveMatrix(1);
		end
	else
		% We're in matlab: the native matlab 2D array can be passed to java as is:
		
		jIntMatrix = octaveMatrix;
	end


end

