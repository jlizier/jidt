% function jDoubleMatrix = octaveToJavaDoubleMatrix(octaveMatrix)
%
% Convert a native octave/matlab matrix to a java double 2D array
%

function jDoubleMatrix = octaveToJavaDoubleMatrix(octaveMatrix)

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		if ((rows(octaveMatrix)*columns(octaveMatrix)) > 1)
			% Do this the normal way
			tmp = javaObject('org.octave.Matrix',reshape(octaveMatrix,1,rows(octaveMatrix)*columns(octaveMatrix)),[rows(octaveMatrix), columns(octaveMatrix)]);
			jDoubleMatrix = tmp.asDoubleMatrix();
		else
			% For length 1 arrays, we need to perform a hack here or else
			%  java thinks the length-one array is a scalar.
			% See octaveToJavaDoubleArray for a further description.
			% So instead we'll do this the slow way (doesn't matter for one element only)
			jDoubleMatrix = javaArray('java.lang.Double', 1, 1);
			jDoubleMatrix(1, 1) = octaveMatrix(1);
		end
	else
		% We're in matlab:
		
		% Presumably there's a quick way to do this in matlab, but since I'm not on matlab, I don't know ...
		% If someone knows or has tested something, please tell me and I'll include it here.
		% In the meantime, we copy element by element
		
		jDoubleMatrix = javaArray('java.lang.Double', rows(octaveMatrix), columns(octaveMatrix));
		for r = 1:rows(octaveMatrix)
			% Slow but effective way:
			for c = 1:columns(octaveMatrix)
				jDoubleMatrix(r,c) = octaveMatrix(r,c);
			end
			% Fast way that doesn't actually work:
			% jDoubleMatrix(r,:) = octaveMatrix(r,:);
		end
	end

end

