% function [calculatedEntropy, winRate, lossRate] = computeEntropyForPlayer(name)
%
% Compute the entropy of moves for a given player, across all games/iterations
% 
% Input:
% - name of the player
%
% Copyright (C) 2017-, Joseph T. Lizier
% Distributed under GNU General Public License v3

function [calculatedEntropy, winRate, lossRate, numGames] = computeEntropyForPlayer(name)

	% Step 1: load all of the player's games' data:
	games = loadGamesForPlayer(name);
	
	% Step 2: the player's moves are in the first column, pull these from
	%  each game into an array of samples that we can compute entropy on:
	moves = [];
	results = [];
	for gameIndex = 1:length(games)
		% Load data from game gameIndex into the variable game
		game = games{gameIndex};
		% First column of game is the player's move, second is opponent's
		%  and third is the result.
		% Pull out the player's moves in this game (first column of game):
		movesInThisGame = ????;
		% Pull out the results in this game (third column of game):
		resultsInThisGame = ????;
		% Append this player's moves to the array we're storing over all iterations:
		moves = [moves; movesInThisGame];
		% Append this player's results to the array over all iterations:
		results = [results; resultsInThisGame];
	end
	
	% Step 3: compute the entropy for this player's moves using our existing scripts:
	calculatedEntropy = ????;
	
	% Step 4: compute the win and loss rates:
	winRate = sum(results == 1)./length(results);
	lossRate = sum(results == -1)./length(results);
	numGames = length(results);
	
	if (nargout == 0)
		fprintf('Entropy for %s over %d iterations: %.4f bits\n', ...
			name, length(moves), calculatedEntropy);
	end
end

