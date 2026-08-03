% function conditionalentropyempirical(xn,yn)
%
% Computes the conditional Shannon entropy over all samples xn of a random
%  variable X, given samples yn of a random variable Y.
%
% Inputs:
% - xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate X.
% - yn - matrix of samples of outcomes x. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate Y.
%        Must have the same number of rows as X.
%
% Outputs:
% - result - conditional Shannon entropy of X given Y
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = conditionalentropyempirical(xn,yn)

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

	% We need to compute H(X,Y) - H(X):
	% 1. joint entropy:
	H_XY = ???;
	% 2. marginal entropy of Y:
	H_Y = ???;
	
	result = H_XY - H_Y;
	
end

