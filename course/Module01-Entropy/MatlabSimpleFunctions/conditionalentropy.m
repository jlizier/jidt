% function conditionalentropy(p)
%
% Computes the conditional Shannon entropy over all outcomes x of a random
%  variable X, given outcomes y of a random variable Y.
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
% - result - conditional Shannon entropy of X given Y
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = conditionalentropy(p)

	% Should we check any potential error conditions on the input?
	% a. Should we check p is a matrix, not a vector?
	% assert(~isvector(p));
	% Actually we won't since a vector would be valid if one variable only ever took one value.
	% b. Check that the probabilities normalise to 1:
	% assert(sum(p(:)) == 1);
	assert(abs(sum(p(:)) - 1) < 0.0000001); % Will work for any dimensionality, and handles numerical rounding errors

	% We need to compute H(X,Y) - H(Y):
	% 1. joint entropy: Can we re-use existing code?
	H_XY = ???;
	% 2. marginal entropy of Y: Can we re-use existing code?
	%  But how to get p_y???
	p_y = ???;
	H_Y = ???;
	
	result = H_XY - H_Y;
	
end

