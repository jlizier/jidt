% function octaveMatrix = javaMatrixToOctave(javaMatrix)
%
% Convert a java matrix (1 or 2D, double or int - but not Integer!!) to an octave or matlab matrix
%
% Octave-java doesn't seem to handle the conversion natively,
%  so we either use org.octave.Matrix (built-in) to do it, or
%  and we must convert each individual array item ourselves (This can be very slow for large matrices)
% 

function octaveMatrix = javaMatrixToOctave(javaMatrix, startRow, startCol, numRows, numCols)

	if (nargin < 2)
		startRow = 1;
	end
	if (nargin < 3)
		startCol = 1;
	end
	if (nargin < 4)
		numRows = size(javaMatrix, 1);
	end
	if (nargin < 5)
		numCols = size(javaMatrix, 2);
	end

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell

		% Convert whole matrix first:
		tmp = javaObject('org.octave.Matrix', javaMatrix);
		% Make sure tmp.ident() is converted to native octave:
		oldFlag = java_convert_matrix (1);
		converted = false;
		unwind_protect
			octaveMatrix = tmp.ident(tmp);
			converted = true;
		unwind_protect_cleanup
			% restore to non-default conversion, otherwise we get
			%  bad errors on other calls
			java_convert_matrix(oldFlag);
		end_unwind_protect
		if (converted)
			if (nargin >= 2)
				% Do some resizing:
				octaveMatrix = octaveMatrix(startRow:startRow+numRows-1, startCol:startCol+numCols-1);
			end
			return;
		end
	else
		% Else we're in matlab, in which case the native java type can be handled, so return it directly:
		octaveMatrix = javaMatrix;
	end

	% Else, we encountered an error in the octave resizing, so fall through to element by element conversion:
	
	octaveMatrix = zeros(numRows, numCols);
	for r = startRow:startRow+numRows-1
		for c = startCol:startCol+numCols-1
			octaveMatrix(r-startRow+1,c-startCol+1) = javaMatrix(r,c);
		end
	end

end

