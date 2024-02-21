# Scissors-Paper-Rock analysis

This set of files are used to analyse the Scissors-Paper-Rock data set.

These are described in full in the tutorial page for this task.

In brief, the files include:
* `setup.m` is called to initialise various system paths for the analysis. You need not call this yourself, the parsing scripts do it for you.
* `listPlayers.m` to pull out the names of players.
* `loadGamesForPlayer.m` to pull out the game data for a given player.
* `translateMove.m` and `translateResult.m` translate the encoded moves and results to text strings.
* `computeEntropyForPlayer.m`, `computeConditionalEntropyForPlayer.m` and `computeMutualInformationForPlayer.m` compute information theoretic measures for a single named player.
* `computeEntropyForAllPlayers.m`, `computeConditionalEntropyForAllPlayers.m` and `computeMutualInformationForAllPlayers.m` compute information-theoretic measures for all players.

