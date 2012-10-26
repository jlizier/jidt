% function jDoubleArray = octaveToJavaDoubleArray(octaveArray)
%
% Convert a native octave array to a java double 1D array
%

function jDoubleArray = octaveToJavaDoubleArray(octaveArray)

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		tmp = javaObject('org.octave.Matrix',octaveArray,[1, length(octaveArray)]);
		jDoubleArray = tmp.asDoubleVector();
	else
		% We're in matlab:
		
		% Presumably there's a quick way to do this in matlab, but since I'm not on matlab, I don't know ...
		% If someone knows or has tested something, please tell me and I'll include it here.
		% In the meantime, we copy element by element
		
		jDoubleArray = javaArray('java.lang.Double', length(octaveArray));
		for r = 1:length(octaveMatrix)
			jDoubleArray(r) = octaveArray(r);
		end
	end

end

