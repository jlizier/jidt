% function infocontent(p)
%
% Computes the Shannon information content for an outcome x of a random variable
%  X with probability p.
%
% Inputs:
% - p - probability to compute the Shannon info content for
%
% Outputs:
% - result - Shannon info content of the probability p
% 
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%

function result = infocontent(p)

	result = log2(1./p);
	
end

