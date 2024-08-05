% function [calculatedMI, winRate, lossRate, numGames] = computeMutualInformationForPlayer(name, fromSelf)
%
% Compute the mutual information of moves for a given player with their own
%  previous move, or the previous move of their opponent
% 
% Input:
% - name (string) name of the player
% - fromSelf (boolean) if true, take MI from the player's own previous move; if false
%    take MI from opponent's previous move.
%
% Copyright (C) 2017-, Joseph T. Lizier
% Distributed under GNU General Public License v3

function [calculatedMI, winRate, lossRate, numGames] = computeMutualInformationForPlayer(name, fromSelf)

	if (nargin < 2)
		fprintf('Defaulting to examine MI from own players past\n');
		fromSelf = true;
	end

	% Step 1: load all of the player's games' data:
	games = loadGamesForPlayer(name);
	
	% Step 2: the player's moves are in the first column, oppenent's in 2nd, pull these from
	%  each game into arrays of samples that we can compute mutual info on:
	nextMoves = [];
	previousMoves = [];
	results = [];
	for gameIndex = 1:length(games)
		% Load data from game gameIndex into the variable game
		game = games{gameIndex};
		% First column of game is the player's move, second is opponent's
		%  and third is the result.
		moves = game(:,1);
		opponentMoves = game(:,2);
		playersResults = game(:,3);
		% Append this player's moves to the array we're storing over all iterations.
		%  TAKE CARE: Can we take all samples here, or only a limited number that
		%  we're able to match up properly to compute mutual information?
		nextMoves = [nextMoves; moves(2:end)];
		if (fromSelf)
			% Grab the previous moves from this player:
			% HINT: This will be the same thing you did in computeConditionalEntropyForPlayer:
			previousMoves = [previousMoves; moves(1:end-1)];
		else
			% Grab the previous moves from their opponent:
			previousMoves = [previousMoves; opponentMoves(1:end-1)];
		end
		% Append this player's results to the array over all iterations:
		%  Which results do we want here -- those of the previous iteration or this one?
		% HINT: This will be the same thing you did in computeConditionalEntropyForPlayer:
		results = [results; playersResults(2:end)];
	end
	
	% Step 3: compute the mutual information for this player's moves using our existing scripts:
	calculatedMI = mutualinformationempirical(nextMoves, previousMoves);
	
	% Step 4: compute the win and loss rates:
	winRate = sum(results == 1)./length(results);
	lossRate = sum(results == -1)./length(results);
	numGames = length(results);
	
	if (nargout == 0)
		fprintf('MI for %s over %d iterations: %.4f bits\n', ...
			name, length(nextMoves), calculatedMI);
	end
end

