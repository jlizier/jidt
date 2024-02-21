% function entropyempirical(xn)
%
% Computes the Shannon entropy over all outcomes x of a random variable
%  X from samples x_n.
%
% Inputs:
% - xn - samples of outcomes x.
%       xn is a column vector, e.g. xn = [0;0;1;0;1;0;1;1;1;0] for a binary variable.
%
% Outputs:
% - result - Shannon entropy over all outcomes
% - symbols - list of unique samples
% - probabilities - probabilities for each sample
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function [result, symbols, probabilities] = entropyempirical(xn)

	% Should we check any potential error conditions on the input?
	assert(isvector(xn));
	%  e.g. what if it is a row vector - can we handle that or
	%   flag error condition? It will work ok!

	% We need to work out the alphabet here.
	% The following returns a vector of the alphabet:
	symbols = unique(xn);
	
	% Next we need to count the number of occurances of each symbol in 
	% the alphabet:
	counts = zeros(1,length(symbols));
	for symbolIndex = 1:length(symbols)
		symbol = symbols(symbolIndex);
		% Count the number of occurances of symbol in xn:
		counts(symbolIndex) = sum(xn == symbol);
	end
	% Now normalise the counts into probabilities:
	probabilities = counts ./ length(xn);
	
	% Once we have probabilities we can simply call our existing function:
	result = entropy(probabilities);
	
end

