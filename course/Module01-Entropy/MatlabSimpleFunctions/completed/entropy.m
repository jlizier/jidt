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
	% assert(sum(p(:)) == 1);
	assert(abs(sum(p(:)) - 1) < 0.0001); % Will work for any dimensionality, and handles numerical rounding errors
	assert(~any(p(:) > 1));
	assert(~any(p(:) < 0));

	% We need to take the expectation value over the Shannon info content at
	%  p(x) for each outcome x:
	
	% Naive: 
	% result = sum(p .* infocontent(p));
	% BUT -- are there any potential error conditions here?
	%  Yes -- if one or more of the values in p is 0!
	
	% Nuanced: (could do for loops here, but will not work if we don't
	%  have the dimensions of p matching the loops).
	% Do p log p first:
	weightedShannonInfos = p .* infocontent(p);
	% Then pick out the p log p values which are not nan
	contributions = weightedShannonInfos(~isnan(weightedShannonInfos));
	% And sum all of them up, over all dimensions:
	result = sum(contributions(:));
	
end

