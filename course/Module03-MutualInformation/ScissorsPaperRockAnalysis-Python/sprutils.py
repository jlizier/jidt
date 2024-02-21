# Basic SPR data plumbing functions
# 
# The following are the basic functions associated to SPR analysis.
# These are based on the original Matlab functions.
#
# Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
# Distributed under GNU General Public License v3

import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import re # This is to easily manipulate the strings in the file name

# Define the global variable for our data path:
dataPath = ""

"""function setDataPath()
Set the data path for where the game files are stored.

Outputs:
- newDataPath - new path to set for the data files

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def setDataPath(newDataPath):
    # uses global variable dataPath
    global dataPath
    dataPath = newDataPath

"""function listPlayers()
Returns a list of all of the players in the scissors-paper-rock data set.

Inputs:
- verbose (default False) - if True print the player names out

Outputs:
- plist - list of all player names

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def listPlayers(verbose=False):
    # uses global variable dataPath
    global dataPath
    
    plist = []
    index = 1
    files = [file for file in listdir(dataPath) if isfile(join(dataPath, file)) & file.endswith('.txt')]
    files.sort()
        
    for file in files:
        # Parse the file name for the player names:
        # a. Pull off the timestamp
        players = [str(x) for x in filter(None, re.split('[,\_,\.]',file))]
        # b. Pull out the players names        
        player1 = players[1]
        player2 = players[2]
        # c. Add them to our list so far
        plist.append(player1)
        plist.append(player2)
        index = index + 2

    # Finally remove any duplicate names:
    plist = list(set(plist))
    plist.sort()
    
    if verbose:
        print('Player names:');
        for name in plist:
            print(name);
    
    return plist

"""function loadGamesForPlayer(name)

Returns a cell array of game sets for the given player name.
For each game set, first column is the player's move, second is their
opponents, and third column is whether they won (1), lost (-1) or drew (0)
 
Inputs:
- name - player name, as a string. Can be '*' to get games for all players
- verbose (default False) - if True print the results are printed to the
 standard output.

Outputs:
- allGameData - list of all games played by this player. Each list item, allGameData[i],
  holds data for a separate game. allGameData[i] is a 2D numpy array where each row
  represents a single iterations within the game. The first column are
  the player's moves (0 == scissors, 1 == paper, 2 == rock), the second
  colummn are the opponents moves, and the third column is the result
  (1 == this player won, 0 == tie, -1 == opponent won).

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def loadGamesForPlayer(name, verbose=False):
    
    # uses global variable dataPath
    global dataPath

    index = 0
    files = [file for file in listdir(dataPath) if isfile(join(dataPath, file)) & file.endswith('.txt')]
    files.sort()
    allGameData = []

    for file in files:
        # Parse the file name for the player names:
        # a. Pull off the timestamp
        players = [str(x) for x in filter(None, re.split('[,\_,\.]',file))]
        # b. Pull out the players names        
        player1 = players[1]
        player2 = players[2]

        if player1 == name or name == '*':
            # Player1 is our player, or we're getting all games
            playerCol = 0
            opponentCol = 1
            thisPlayer = player1
            opponent = player2      
        elif player2 == name:
            # Player2 is our player
            playerCol = 1
            opponentCol = 0
            thisPlayer = player2
            opponent = player1
        else:
            continue # Move to next file

        # Load this game in:
        df = pd.read_csv(join(dataPath, file),sep="\t",comment="%", header=None)
        gameData = df.values
        # Grab their moves:
        # 0 = scissors
        # 1 = paper
        # 2 = rock
        playerMoves = gameData[:,playerCol]
        opponentMoves = gameData[:,opponentCol]
		# The player wins if their move is one
		#  less than opponents, or (opponent - player) mod 3 == 1.
		# If (opponent - player) mod 3 == 2, then opponent wins.
		# Otherwise if (player == opponent) then it's a tie.
		# Can express this concisely as the following to make
		#  I win == 1
		#  You win == -1
		#  Tie == 0
        results = ((opponentMoves - playerMoves + 1) % 3) - 1
		# Now store all of this as a new entry in the list:
        allGameData.append(np.column_stack((playerMoves, opponentMoves, results)))

        if verbose:
            # User wants the games printed:
            print("Game {} for {} ({} iterations):".format(index, name, allGameData[index].shape[0]))
            for i in range(allGameData[index].shape[0]):
                # allGameData[index][i,:] is the data for this one iteration
                # in the game
                print('{}:\t{},\t{}:\t{},\tresult: {}'.format(thisPlayer, translateMove(allGameData[index][i,0]), opponent, translateMove(allGameData[index][i,1]), translateResult(allGameData[index][i,2])))
            print()

        if name == '*':
            # If we're grabbing data for all players, then take the 
            #  player2's perspective as well:
            index = index + 1
            allGameData.append(np.column_stack((opponentMoves, playerMoves, -results)))
        
            if verbose:
                # User wants the games printed:
                print("Game {} for {} ({} iterations):".format(index, name, allGameData[index].shape[0]))
                for i in range(allGameData[index].shape[0]):
                    # allGameData[i,:] is the data for this one iteration
                    # in the game
                    print('{}:\t{},\t{}:\t{},\tresult: {}'.format(opponent, translateMove(allGameData[index][i,0]), thisPlayer, translateMove(allGameData[index][i,1]), translateResult(allGameData[index][i,2])))

        index += 1

    if (index == 0): raiseError('No games found for user {}'.format(name))
    
    return allGameData


"""function translateMove(move)

Returns a string representation of the given move index:
0 -> scissors
1 -> paper
2 -> rock

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def translateMove(move):
    if move == 0:
        stringRepresentation = 'scis'
    elif move == 1:
        stringRepresentation = 'papr'
    elif move == 2:
        stringRepresentation = 'rock'
    else:
        print('Error: Move {} not recognised'.format(move))

    return stringRepresentation

"""function translateResult(gameResult)

Returns a string representation of the given game result:
-1 -> lose
0 -> tie
1 -> win

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def translateResult(gameResult):
    if gameResult == -1:
        stringRepresentation = 'los'
    elif gameResult == 0:
        stringRepresentation = 'tie'
    elif gameResult == 1:
        stringRepresentation = 'win'
    else:
        print('Error: Result {} not recognised'.format(gameResult))

    return stringRepresentation

