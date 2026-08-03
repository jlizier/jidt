# Scissors-Paper-Rock analysis

This set of files are used to analyse the Scissors-Paper-Rock data set in this tutorial task.

The Matlab live script `ScissorsPaperRockAnalysis.mlx` guides you through the task and describes the role of each script which it uses.

In brief, these other scripts include:
* `listPlayers.m` to pull out the names of players.
* `loadGamesForPlayer.m` to pull out the game data for a given player.
* `translateMove.m` and `translateResult.m` translate the encoded moves and results to text strings.
* `computeEntropyForPlayer.m`, `computeConditionalEntropyForPlayer.m` and `computeMutualInformationForPlayer.m` compute information theoretic measures for a single named player -- these will be completed during the tutorial task.
* `computeEntropyForAllPlayers.m`, `computeConditionalEntropyForAllPlayers.m` and `computeMutualInformationForAllPlayers.m` compute information-theoretic measures for all players -- these will be completed during the tutorial task.
* `setup.m` is not used by the live script, but is left here for uses wishing to analyse the data outside of the live script. It can be called to initialise various system paths for the analysis before calling the main scripts.

