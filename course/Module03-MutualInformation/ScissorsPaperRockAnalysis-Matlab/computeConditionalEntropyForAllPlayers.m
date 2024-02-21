% function [names, entropies, winRates, lossRates] = computeConditionalEntropyForAllPlayers()
%
% Compute the conditional entropy of moves for each player, conditioned on their previous move,
%  across all games/iterations
% 

function [names, entropies, winRates, lossRates] = computeConditionalEntropyForAllPlayers()

	% Step 1: load all of the player's names:
	names = listPlayers();

	% Step 2: compute conditional entropy for each player:
	index = 1;
	entropies = zeros(length(names),1);
	winRates = zeros(length(names),1);
	lossRates = zeros(length(names),1);
	for name = names
		% Compute the entropy for the moves of this player.
		% HINT: use the script that you just completed; the player's 
		%  name as a string to pass in is name{:}.
		[calculatedEntropy, winRate, lossRate, numGames] = ...
			???;
		fprintf('%s: %.4f bits,\twin rate = %.4f,\tloss rate = %.4f, num games = %d\n', ...
			name{:}, calculatedEntropy, winRate, lossRate, numGames);
		
		entropies(index) = calculatedEntropy;
		winRates(index) = winRate;
		lossRates(index) = lossRate;

		index = index + 1;
	end

	% Plot the winRates and lossRates versus entropies:
	figure(1);
	plot(entropies, winRates, 'x');
	title('Win rates versus cond entropies of single players');
	xlabel('Entropy of moves (bits)');
	ylabel('Win rate');
	figure(2);
	plot(entropies, lossRates, 'x');
	title('Loss rates versus cond entropies of single players');
	xlabel('Entropy of moves (bits)');
	ylabel('Loss rate');
	
	% Compute correlations of entropy to win rate and to loss rate:
	winToEntropyCorr = ???;
	lossToEntropyCorr = ???;
	fprintf('Correlation of win  rate to cond entropy is: %.4f\n', ...
		winToEntropyCorr);
	fprintf('Correlation of loss rate to cond entropy is: %.4f\n', ...
		lossToEntropyCorr);
	% Are these statistically significant?
	% Can you adjust your code to check for that?
	% HINT: look at other return values from the correlation function
	
end

