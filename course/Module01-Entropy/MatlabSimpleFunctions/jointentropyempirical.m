% function jointentropyempirical(xn, yn)
%
% Computes the Shannon entropy over all outcome vectors x of a vector random
%  variable X from sample vectors x_n. User can call with two such arguments 
%  if they don't wish to join them outside of the call.
%
% Inputs:
% - xn - matrix of samples of outcomes x. May be a 1D vector of samples
%        (in which case yn is also supplied), or
%        a 2D matrix, where each row is a vector sample for a multivariate X
%        (in which case yn is not supplied).
%  - yn - as per xn, except that yn is not required to be supplied (in which
%	 case the entropy is only calculated over the multivariate xn variable).
%
% Outputs:
% - result - joint Shannon entropy over all samples
% - symbols - list of unique joint vector samples
% - probabilities - probabilities for each joint symbol
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function [result, symbols, probabilities] = jointentropyempirical(xn, yn)

	% Should we check any potential error conditions on the input?
	assert(length(size(xn))==2);
	% Convert to column vectors if not already:
	if (size(xn,1) == 1)
		% xn has only one row, assume these are multiple observations of single dimensional variable:
		xn = xn'; % Transpose it so it is only column
	end
	if (nargin > 1)
		% Two arguments
		assert(length(size(yn))==2);
		if (size(yn,1) == 1)
			% yn has only one row, assume these are multiple observations of single dimensional variable
			yn = yn'; % Transpose it so it is only column
		end
		% Check that their number of rows are the same:
		assert(size(xn,1) == size(yn,1));
		% Now joint them up so we only need work with xn
		xn = [xn,yn]; % Joins the column vectors into a matrix
	end
	% Now, we are only working with a 2D matrix xn of row vector samples
	%   (i.e. each column represents a variable/dimension, while each row is
	%    one sample of the joint variable)
	
	% TRICK: Next combine the row vectors in each sample into a single 
	%  symbol, so that we can simply compute entropy on that combined symbol
	[symbols,~,combinedSamples] = unique(xn, 'rows');
	
	% And return the entropy:
	[result, ~, probabilities] = entropyempirical(combinedSamples);

	% The order of symbols is the same as their order for the probabilities
    
end

