% function mutualinformation(p)
%
% Computes the mutual information over all outcomes x of a random
%  variable X with outcomes y of a random variable Y.
%  Probability matrix p(x,y) is given for each candidate outcome
%  (x,y).
%
% Inputs:
% - p - 2D probability distribution function over all outcomes (x,y).
%       p is a matrix over all combinations of x and y,
%	where p(1,3) gives the probability of the first symbol of variable
%	x co-occuring with the third symbol of variable y.
%       E.g. p = [0.2, 0.3; 0.1, 0.4]. The sum over p must be 1.
%
% Outputs:
% - result - mutual information of X with Y
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = mutualinformation(p)

	% Should we check any potential error conditions on the input?
	% a. Should we check p is a matrix, not a vector?
	% assert(~isvector(p));
	% Actually we won't since a vector would be valid if one variable only ever took one value
	% b. Check that the probabilities normalise to 1:
	% assert(sum(p(:)) == 1);
	assert(abs(sum(p(:)) - 1) < 0.0000001); % Will work for any dimensionality, and handles numerical rounding errors

	% We need to compute H(X) + H(Y) - H(X,Y):
	% 1. joint entropy:
	H_XY = jointentropy(p);
	% 2. marginal entropy of X:
	%  But how to get p_x???
	p_x = sum(p,2); % Since x changes along the rows, summing over the y's (dimension 2 argument in the sum) will just return p(x)
	H_X = entropy(p_x);
	% 2. marginal entropy of Y:
	%  But how to get p_y???
	p_y = sum(p,1); % Since y changes along the columns, summing over the x's (dimension 1 argument in the sum) will just return p(y)
	H_Y = entropy(p_y);
	
	result = H_X + H_Y - H_XY;
	
end

