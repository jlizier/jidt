% function jDoubleArray = octaveToJavaDoubleArray(octaveArray)
%
% Convert a native octave array to a java double 1D array
%

function jDoubleArray = octaveToJavaDoubleArray(octaveArray)

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		if (length(octaveArray) > 1)
			% Do this the normal way
			tmp = javaObject('org.octave.Matrix',octaveArray,[1, length(octaveArray)]);
			jDoubleArray = tmp.asDoubleVector();
		else
			% For length 1 arrays, we need to perform a hack here or else
			%  java thinks the length-one array is a scalar.
			% I thought I had this work once:
			%  tmp = javaObject('org.octave.Matrix',[octaveArray,octaveArray] ,[1, 1]);
			%  jDoubleArray = tmp.asDoubleVector();
			% but now can't get that to repeat.
			% So instead we'll do this the slow way (doesn't matter for one element only)
			jDoubleArray = javaArray('java.lang.Double', 1);
			jDoubleArray(1) = octaveArray(1);
		end
	else
		% We're in matlab: the native matlab array can be passed to java as is:
		
		jDoubleArray = octaveArray;
	end

end

