% function mutualinformationempirical(xn,yn)
%
% Computes the mutual information over all samples xn of a random
%  variable X with samples yn of a random variable Y.
%
% Inputs:
% - xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate X.
% - yn - matrix of samples of outcomes x. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate Y.
%        Must have the same number of rows as X.
%
% Outputs:
% - result - mutual information of X with Y
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = mutualinformationempirical(xn,yn)

	% Should we check any potential error conditions on the input?
	if (isvector(xn))
		% Convert it to column vector if not already:
		if (size(xn,1) == 1)
			% xn has only one row:
			xn = xn'; % Transpose it so it is only column
		end
	end
	if (isvector(yn))
		% Convert it to column vector if not already:
		if (size(yn,1) == 1)
			% yn has only one row:
			yn = yn'; % Transpose it so it is only column
		end
	end
	% Check that their number of rows are the same:
	assert(size(xn,1) == size(yn,1));

	% We need to compute H(X) + H(Y) - H(X,Y):
	% 1. joint entropy:
	H_XY = ???; % How to compute this empirically ...?
	% 2. marginal entropy of Y: (calling 'joint' in case yn is multivariate)
	H_Y = ???;
	% 3. marginal entropy of X: (calling 'joint' in case yn is multivariate)
	H_X = ???;
	
	result = H_X + H_Y - H_XY;
end

