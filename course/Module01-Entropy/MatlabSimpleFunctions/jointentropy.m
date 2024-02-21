% function jointentropy(p)
%
% Computes the joint Shannon entropy over all outcome vectors x of a vector
%  random variable X with probability matrix p(x) for each candidate outcome
%  vector x.
%
% Inputs:
% - p - probability distribution function over all outcome vectors x.
%       p is a matrix over all combinations of the sub-variables of x,
%	where p(1,3) gives the probability of the first symbol of sub-variable
%	x1 co-occuring with the third symbol of sub-variable x2.
%       E.g. p = [0.2, 0.3; 0.1, 0.4]. The sum over p must be 1.
%
% Outputs:
% - result - joint Shannon entropy of the probability distribution p
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = jointentropy(p)

	% Should we check any potential error conditions on the input?

	% We need to take the expectation value over the Shannon info content at
	%  p(x) for each outcome x in the joint PDF:
	% Hint: will your code for entropy(p) work, or can you alter it slightly
	%  to make it work?
	???
	
end

