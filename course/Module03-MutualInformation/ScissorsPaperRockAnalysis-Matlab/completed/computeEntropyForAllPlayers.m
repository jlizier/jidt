% function [names, entropies, winRates, lossRates] = computeEntropyForAllPlayers()
%
% Compute the entropy of moves for each player, across all games/iterations
% 

function [names, entropies, winRates, lossRates] = computeEntropyForAllPlayers()

	% Step 1: load all of the player's names:
	names = listPlayers();

	% Step 2: compute entropy for each player:
	index = 1;
	entropies = zeros(length(names),1);
	winRates = zeros(length(names),1);
	lossRates = zeros(length(names),1);
	for name = names
		[calculatedEntropy, winRate, lossRate, numGames] = ...
			computeEntropyForPlayer(name{:});
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
	title('Win rates versus entropies of single players');
	xlabel('Entropy of moves (bits)');
	ylabel('Win rate');
	figure(2);
	plot(entropies, lossRates, 'x');
	title('Loss rates versus entropies of single players');
	xlabel('Entropy of moves (bits)');
	ylabel('Loss rate');
	
	% Compute correlations and check if these are statistically significant:
	% Are these statistically significant?
	if (exist ('OCTAVE_VERSION', 'builtin'))
		% This is running on Octave (not Matlab), so do this the hard way:
		winToEntropyCorr = corr(winRates, entropies);
		lossToEntropyCorr = corr(lossRates, entropies);
		% Now compute the pValues:
		winToEntropyCorrTValue = winToEntropyCorr ./ ...
					sqrt((1-winToEntropyCorr.^2) ./ (length(names)-2));
		lossToEntropyCorrTValue = lossToEntropyCorr ./ ...
					sqrt((1-lossToEntropyCorr.^2) ./ (length(names)-2));
		% Using two-tailed tests:
		winToEntropyCorrTCdf = tcdf(winToEntropyCorrTValue, length(names)-2);
		if (winToEntropyCorrTCdf < 0.5)
			% Account for the probability mass on the other tail of the distribution:
			winToEntropyCorrPValue = winToEntropyCorrTCdf .* 2;
		else
			% Account for the probability mass on the other tail of the distribution:
			winToEntropyCorrPValue = 2.*(1 - winToEntropyCorrTCdf);
		end
		lossToEntropyCorrTCdf = tcdf(lossToEntropyCorrTValue, length(names)-2);
		if (lossToEntropyCorrTCdf < 0.5)
			% Account for the probability mass on the other tail of the distribution:
			lossToEntropyCorrPValue = lossToEntropyCorrTCdf .* 2;
		else
			% Account for the probability mass on the other tail of the distribution:
			lossToEntropyCorrPValue = 2.*(1 - lossToEntropyCorrTCdf);
		end
	else
		% We're running on Matlab, so do this the easy way:
		[winToEntropyCorr, winToEntropyCorrPValue] = corr(winRates, entropies);
		[lossToEntropyCorr, lossToEntropyCorrPValue] = corr(lossRates, entropies);
	end
	fprintf('Correlation of win  rate to entropy is: %.4f (pValue %.4f)\n', ...
		winToEntropyCorr, winToEntropyCorrPValue);
	fprintf('Correlation of loss rate to entropy is: %.4f (pValue %.4f)\n', ...
		lossToEntropyCorr, lossToEntropyCorrPValue);
	
end

