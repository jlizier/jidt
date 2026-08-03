% function [names, mutualInfos, winRates, lossRates] = computeMutualInformationForAllPlayers(fromSelf)
%
% Compute the mutual information of moves for each player with their own
%  previous move, or the previous move of their opponent, across all games/iterations.
% 
% Input:
% - fromSelf (boolean) if true, take MI from the player's own previous move; if false
%    take MI from opponent's previous move.
%
% Copyright (C) 2017-, Joseph T. Lizier
% Distributed under GNU General Public License v3

function [names, mutualInfos, winRates, lossRates] = computeMutualInformationForAllPlayers(fromSelf)

	if (nargin < 1)
		fromSelf = true;
	end

	% Step 1: load all of the player's names:
	names = listPlayers();

	% Step 2: compute mutual info for each player:
	index = 1;
	mutualInfos = zeros(length(names),1);
	winRates = zeros(length(names),1);
	lossRates = zeros(length(names),1);
	for name = names
        % Compute the mutual info for the moves of this player.
        % HINT: use the script that you just completed passing in name and fromSelf
		[calculatedMI, winRate, lossRate, numGames] = ...
			computeMutualInformationForPlayer(name{:}, fromSelf);
		fprintf('%s: %.4f bits,\twin rate = %.4f,\tloss rate = %.4f, num games = %d\n', ...
			name{:}, calculatedMI, winRate, lossRate, numGames);
		
		mutualInfos(index) = calculatedMI;
		winRates(index) = winRate;
		lossRates(index) = lossRate;

		index = index + 1;
	end

	% Plot the winRates and lossRates versus mutualInfos:
	figure(1);
	plot(mutualInfos, winRates, 'x');
	title('Win rates versus mutual information for single players');
	xlabel('Mutual information of moves (bits)');
	ylabel('Win rate');
	figure(2);
	plot(mutualInfos, lossRates, 'x');
	title('Loss rates versus mutual information for single players');
	xlabel('Mutual information of moves (bits)');
	ylabel('Loss rate');
	
	% Compute correlations and check if these are statistically significant:
	% Are these statistically significant?
	if (exist ('OCTAVE_VERSION', 'builtin'))
		% This is running on Octave (not Matlab), so do this the hard way:
		winToMICorr = corr(winRates, mutualInfos);
		lossToMICorr = corr(lossRates, mutualInfos);
		% Now compute the pValues:
		winToMICorrTValue = winToMICorr ./ ...
					sqrt((1-winToMICorr.^2) ./ (length(names)-2));
		lossToMICorrTValue = lossToMICorr ./ ...
					sqrt((1-lossToMICorr.^2) ./ (length(names)-2));
		% Using two-tailed tests:
		winToMICorrTCdf = tcdf(winToMICorrTValue, length(names)-2);
		if (winToMICorrTCdf < 0.5)
			% Account for the probability mass on the other tail of the distribution:
			winToMICorrPValue = winToMICorrTCdf .* 2;
		else
			% Account for the probability mass on the other tail of the distribution:
			winToMICorrPValue = 2.*(1 - winToMICorrTCdf);
		end
		lossToMICorrTCdf = tcdf(lossToMICorrTValue, length(names)-2);
		if (lossToMICorrTCdf < 0.5)
			% Account for the probability mass on the other tail of the distribution:
			lossToMICorrPValue = lossToMICorrTCdf .* 2;
		else
			% Account for the probability mass on the other tail of the distribution:
			lossToMICorrPValue = 2.*(1 - lossToMICorrTCdf);
		end
	else
		% We're running on Matlab, so do this the easy way:
		[winToMICorr, winToMICorrPValue] = corr(winRates, mutualInfos);
		[lossToMICorr, lossToMICorrPValue] = corr(lossRates, mutualInfos);
	end
	fprintf('Correlation of win  rate to MI is: %.4f (pValue %.4f)\n', ...
		winToMICorr, winToMICorrPValue);
	fprintf('Correlation of loss rate to MI is: %.4f (pValue %.4f)\n', ...
		lossToMICorr, lossToMICorrPValue);
	
end

