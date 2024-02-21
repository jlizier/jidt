% function games = loadGamesForPlayer(name)
%
% Returns a cell array of game sets for the given player name.
% For each game set, first column is the player's move, second is their
% opponents, and third column is whether they won (1), lost (-1) or drew (0)
% 
% Inputs:
% - name - player name, as a string. Can be '*' to get games for all players
% Outputs:
% - games - cell array of all games played by this player. Each cell, games{i},
%   holds data for a separate game. Each game{i} is a 2D array where each row
%   represents a single iterations within the game. The first column are
%   the player's moves (0 == scissors, 1 == paper, 2 == rock), the second
%   colummn are the opponents moves, and the third column is the result
%   (1 == this player won, 0 == tie, -1 == opponent won).
% If no output argument is requested, the results are simply printed to the
%  standard output.
%  
function games = loadGamesForPlayer(name)

	setup

	files = dir([dataPath, '*.txt']);
	index = 1;
	allGameData = {};
	for file = files'
		% Parse the file name for the player names:
		% a. Pull off the timestamp
		[timestamp, remainder] = strtok(file.name, '_');
		[player1, remainder] = strtok(remainder, '_');
		remainder(1) = []; % Remove the leading '_' (not sure why this isn't required above ...)
		[player2, remainder] = strtok(remainder, '.');
		% fprintf('player1: %s, player2: %s\n', player1, player2);
		
		if (strcmp(player1, name) || strcmp('*', name))
			% Player1 is our player, or we're getting all games
			playerCol = 1;
			opponentCol = 2;
			thisPlayer = player1;
			opponent = player2;
		elseif (strcmp(player2, name))
			% Player2 is our player
			playerCol = 2;
			opponentCol = 1;
			thisPlayer = player2;
			opponent = player1;
		else
			continue; % Move to next file
		end
		% Load this game in:
		gameData = load([dataPath, file.name]);
		% Grab their moves:
		% 0 = scissors
		% 1 = paper
		% 2 = rock
		playerMoves = gameData(:,playerCol);
		opponentMoves = gameData(:,opponentCol);
		% The player wins if their move is one
		%  less than opponents, or (opponent - player) mod 3 == 1.
		% If (opponent - player) mod 3 == 2, then opponent wins.
		% Otherwise if (player == opponent) then it's a tie.
		% Can express this concisely as the following to make
		%  I win == 1
		%  You win == -1
		%  Tie == 0
		results = mod(opponentMoves - playerMoves + 1, 3) - 1;
		% Now store all of this in the cell array:
		allGameData{index} = [playerMoves, opponentMoves, results];

		if (nargout == 0)
			% User doesn't want the data returned, just printed:
			fprintf('Game %d for %s (%d iterations):\n', ...
				index, name, size(allGameData{index},1));
			for iteration = allGameData{index}'
				% iteration is the data for this one iteration
				%  in the game
				fprintf('%s:\t%s,\t%s:\t%s,\tresult: %s\n', ...
					thisPlayer, translateMove(iteration(1)), ...
					opponent, translateMove(iteration(2)), ...
					translateResult(iteration(3)));
			end
		end

		if (strcmp('*', name))
			% If we're grabbing data for all players, then take the 
			%  player2's perspective as well:
			index = index + 1;
			allGameData{index} = [opponentMoves, playerMoves, -results];			
			if (nargout == 0)
				% User doesn't want the data returned, just printed:
				fprintf('Game %d for %s (%d iterations):\n', ...
					index, name, size(allGameData{index},1));
				for iteration = allGameData{index}'
					% iteration is the data for this one iteration
					%  in the game
					fprintf('%s:\t%s,\t%s:\t%s,\tresult: %s\n', ...
						player2, translateMove(iteration(1)), ...
						player1, translateMove(iteration(2)), ...
						translateResult(iteration(3)));
				end
			end
		end
				
		index = index + 1;
	end
	
	if (nargout ~= 0)
		% User wants the game data returned:
		games = allGameData;
	end

	if (index == 1)
		error('No games found for user %s', name);
	end
end


