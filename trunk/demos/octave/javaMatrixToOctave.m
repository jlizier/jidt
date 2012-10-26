% function octaveMatrix = javaMatrixToOctave(javaMatrix)
%
% Convert a java matrix (1 or 2D, double or int - but not Integer!!) to an octave matrix
%
% Unfortunately octave-java doesn't seem to handle the conversion properly,
%  and we must convert each individual array item ourselves (This can be very slow for large matrices)
%  or with org.octave.Matrix (built-in).
% 

function octaveMatrix = javaMatrixToOctave(javaMatrix, startRow, startCol, numRows, numCols)

	if (nargin < 2)
		startRow = 1;
	end
	if (nargin < 3)
		startCol = 1;
	end
	if (nargin < 4)
		numRows = rows(javaMatrix);
	end
	if (nargin < 5)
		numCols = columns(javaMatrix);
	end

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		%  (only use if we're converting the whole array though, too many issues
		%  with type checking etc to bother otherwise)

		if (nargin < 2)
			tmp = javaObject('org.octave.Matrix', javaMatrix);
			% Make sure tmp.ident() is converted to native octave:
			oldFlag = java_convert_matrix (1);
			unwind_protect
				octaveMatrix = tmp.ident(tmp);
			unwind_protect_cleanup
				% restore to non-default conversion, otherwise we get
				%  bad errors on other calls
				java_convert_matrix(oldFlag);
			end_unwind_protect
			return;
		end
	end

	% Else, either we couldn't be bothered doing the resizing for java,
	%  or we were in Matlab all along. (If there's a fast way for matlab, tell me)
	
	octaveMatrix = zeros(numRows, numCols);
	for r = startRow:startRow+numRows-1
		for c = startCol:startCol+numCols-1
			octaveMatrix(r-startRow+1,c-startCol+1) = javaMatrix(r,c);
		end
	end

end

