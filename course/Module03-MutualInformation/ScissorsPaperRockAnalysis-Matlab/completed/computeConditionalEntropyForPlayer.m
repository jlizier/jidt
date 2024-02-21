% function [calculatedEntropy, winRate, lossRate] = computeConditionalEntropyForPlayer(name)
%
% Compute the conditional entropy of moves for a given player, conditioned on
%  their previous move across all games/iterations
% 
% Input:
% - name of the player

function [calculatedEntropy, winRate, lossRate, numGames] = computeConditionalEntropyForPlayer(name)

	% Step 1: load all of the player's games' data:
	games = loadGamesForPlayer(name);
	
	% Step 2: the player's moves are in the first column, pull these from
	%  each game into arrays of samples that we can compute conditional entropy on:
	nextMoves = [];
	previousMoves = [];
	results = [];
	for gameIndex = 1:length(games)
		% Load data from game gameIndex into the variable game
		game = games{gameIndex};
		% First column of game is the player's move, second is opponent's
		%  and third is the result.
		moves = game(:,1);
		playersResults = game(:,3);
		% Append this player's moves to the array we're storing over all iterations.
		%  TAKE CARE: Can we take all samples here, or only a limited number that
		%  we're able to match up properly to compute conditional entropy?
		nextMoves = [nextMoves; moves(2:end)];
		previousMoves = [previousMoves; moves(1:end-1)];
		% Append this player's results to the array over all iterations:
		%  Which results do we want here -- those of the previous iteration or this one?
		results = [results; playersResults(2:end)];
	end
	
	% Step 3: compute the condtional entropy for this player's moves using our existing scripts:
	calculatedEntropy = conditionalentropyempirical(nextMoves, previousMoves);
	
	% Step 4: compute the win and loss rates:
	winRate = sum(results == 1)./length(results);
	lossRate = sum(results == -1)./length(results);
	numGames = length(results);
	
	if (nargout == 0)
		fprintf('Conditional entropy for %s over %d iterations: %.4f bits\n', ...
			name, length(nextMoves), calculatedEntropy);
	end
end

