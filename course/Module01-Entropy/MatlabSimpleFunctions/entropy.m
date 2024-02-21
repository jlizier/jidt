% function entropy(p)
%
% Computes the Shannon entropy over all outcomes x of a random variable
%  X with probability vector p(x) for each candidate outcome x.
%
% Inputs:
% - p - probability distribution function over all outcomes x.
%       p is a vector, e.g. p = [0.25, 0.75], the sum over which must be 1.
%
% Outputs:
% - result - Shannon entropy of the probability distribution p
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = entropy(p)

	% Should we check any potential error conditions on the input?

	% We need to take the expectation value over the Shannon info content at
	%  p(x) for each outcome x:
 	%Alter the equation below to provide the correct entropy:
	result = ???;
	
end

