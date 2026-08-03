% function conditionalmutualinformation(p)
%
% Computes the mutual information over all outcomes x of a random
%  variable X with outcomes y of a random variable Y, conditioning on 
%  outcomes z of a random variable Z.
%  Probability matrix p(x,y,z) is given for each candidate outcome
%  (x,y,z).
%
% Inputs:
% - p - 3D probability distribution function over all outcomes (x,y,z).
%       p is a matrix over all combinations of x and y and z,
%	where p(1,3,2) gives the probability of the first symbol of variable
%	x co-occuring with the third symbol of variable y and the second
%	symbol of z.
%       The sum over p must be 1.
%       E.g.:
%         p(:,:,1) = [0.114286, 0.171429; 0.057143, 0.228571];
%         p(:,:,2) = [0.171429, 0.114286; 0.028571, 0.114286];
%
% Outputs:
% - result - mutual information of X with Y given Z
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = conditionalmutualinformation(p)

	% Should we check any potential error conditions on the input?
	% a. Should we check p is a 3D matrix, not a vector or 2D matrix?
	% assert(~isvector(p));
	% Actually we won't since a vector/2D matrix would be valid if one/two variable only ever took one value
	% b. Check that the probabilities normalise to 1:
	% assert(sum(p(:)) == 1);
	assert(abs(sum(p(:)) - 1) < 0.0000001); % Will work for any dimensionality, and handles numerical rounding errors

	% We need to compute H(X|Z) + H(Y|Z) - H(X,Y|Z).
	% But our conditional entropy calculator won't do H(X,Y|Z) since it doesn't accept a joint probability for X,Y.
	% So, easier to rewrite as:
	%  H(X,Z) - H(Z) + H(Y,Z) - H(Z) - H(X,Y,Z) + H(Z)
	%  = H(X,Z) - H(Z) + H(Y,Z) - H(X,Y,Z)
	% 1. joint entropy:
	H_XYZ = jointentropy(p);
	% 2. entropy of X,Z:
	%  But how to get p_xz???
	% Sum p over the y's (dimension 2 argument in the sum) will just return p(x,z) terms. Won't be a 2D array, but fine to compute entropy on
	p_xz = sum(p,2);
	H_XZ = jointentropy(p_xz);
	% 3. entropy of Y,Z:
	%  But how to get p_yz???
	% Sum p over the x's (dimension 1 argument in the sum) will just return p(y,z) terms. Won't be a 2D array, but fine to compute entropy on
	p_yz = sum(p,1);
	H_YZ = jointentropy(p_yz);
	% 4. marginal entropy of Z:
	%  But how to get p_z???
	% Sum p_xz over the x's (dimension 1 argument in the sum) will just return p(z) terms. Won't be a 1D array, but fine to compute entropy on
	p_z = sum(p_xz,1); 
	H_Z = jointentropy(p_z);
	
	result = H_XZ - H_Z + H_YZ - H_XYZ;
	
end

