% function playerslist = listPlayers()
%
% If called as:
% > listPlayers
% i.e. without a return argument, writes a list of all of the players in the
%  scissors-paper-rock data set.
%
% Otherwise, if called as:
% > playerslist = listPlayers();
% i.e. with a return argument, returns a cell array of all of the players in
% the scissors-paper-rock data set.
% Each element can be accessed from the cell array using e.g. playerslist{i}
% 
% Assumes that dataPath is defined in setup.m
%
% Copyright (C) 2017-, Joseph T. Lizier
% Distributed under GNU General Public License v3

function playerslist = listPlayers()

    % When running from the live script:
	global dataPath
    % When running from command line:
    % setup
    
	plist = {};

	files = dir([dataPath, '*.txt']);
	index = 1;
	for file = files'
		% Parse the file name for the player names:
		% a. Pull off the timestamp
		[timestamp, remainder] = strtok(file.name, '_');
		[player1, remainder] = strtok(remainder, '_');
		remainder(1) = []; % Remove the leading '_' (not sure why this isn't required above ...)
		[player2, remainder] = strtok(remainder, '.');
		% fprintf('player1: %s, player2: %s\n', player1, player2);
		plist{index} = player1;
		plist{index+1} = player2;
		index = index + 2;
	end

	% Finally remove any duplicate names:
	plist = unique(plist);
	
	if (nargout == 0)
		% The user did not ask for a return value here		
		fprintf('Player names:\n');
		for name = plist
			fprintf('%s\n', name{:});
		end
	else
		playerslist = plist;
	end
end

