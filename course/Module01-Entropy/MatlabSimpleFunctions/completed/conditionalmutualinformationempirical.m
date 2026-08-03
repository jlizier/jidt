% function conditionalmutualinformationempirical(xn,yn,zn)
%
% Computes the mutual information over all samples xn of a random
%  variable X with samples yn of a random variable Y, conditioning on 
%  samples zn of a random variable Z.
%
% Inputs:
% - xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate X.
% - yn - matrix of samples of outcomes y. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate Y.
%        Must have the same number of rows as X.
% - zn - matrix of samples of outcomes z. May be a 1D vector of samples, or
%        a 2D matrix, where each row is a vector sample for a multivariate Z
%        which will be conditioned on.
%        Must have the same number of rows as X.
%
% Outputs:
% - result - conditional mutual information of X with Y, given Z
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = conditionalmutualinformationempirical(xn,yn,zn)

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
	if (isvector(zn))
		% Convert it to column vector if not already:
		if (size(zn,1) == 1)
			% zn has only one row:
			zn = zn'; % Transpose it so it is only column
		end
	end
	% Check that their number of rows are the same:
	assert(size(xn,1) == size(yn,1));
	assert(size(xn,1) == size(zn,1));

	% We need to compute H(X|Z) + H(Y|Z) - H(X,Y|Z):
	% 1. conditional joint entropy:
	H_XY_given_Z = conditionalentropyempirical([xn, yn], zn);
	% 2. conditional entropy of Y:
	H_Y_given_Z = conditionalentropyempirical(yn, zn);
	% 3. conditional entropy of X:
	H_X_given_Z = conditionalentropyempirical(xn, zn);
	
	result = H_X_given_Z + H_Y_given_Z - H_XY_given_Z;

	% Alternatively, note that we could compute I(X;Y,Z) - I(X;Z)
	% 1. joint MI:
	% I_X_YZ = mutualinformationempirical(xn, [yn, zn]);
	% 2. MI just from Z:
	% I_X_Z = mutualinformationempirical(xn, zn);
	% Then:
	% result = I_X_YZ - I_X_Z;

end

