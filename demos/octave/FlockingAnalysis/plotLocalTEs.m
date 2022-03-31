function plotLocalTEs(properties)
% Plot the local TEs to show where the information transfer hotspots are from target fish relative to each source
%
% Author: Joseph T. Lizier, 2019
%
% Inputs:
% - properties - object with properties for the calculations,
%   with sub-members as specificied in the loadProperties.m file. If not supplied
%   the properties are loaded from loadProperties.m

%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2019, Joseph T. Lizier et al.
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

	if (nargin < 1)
		fprintf('No properties object supplied, attempting to load properties via a loadProperties script ...');
		% By default, just try to load properties locally
		if (exist('loadProperties') == 2)
			% there is a loadProperties script
			loadProperties;
		else
			% there is not a loadProperties script
			error('No properties object supplied, and no loadProperties script found.');
		end
	end

	load(properties.resultsFile);
	% Loads:
	% S -- source samples
	% D -- target samples
	% Dpast -- target past samples (embedded up to k previous samples)
	% files -- cell array of file names that we took samples from
	% fileTimeAndPair -- each row holds file index, time index, target index, source index
	% RelSourcePos -- each row holds distance between the pair for this sample,
	%	their xyRelativeAngleOfSource, and zRelativeAngleOfSource
	% lag -- source-target lag that is in use
	% k -- embedding length for target that is in use
	% tau -- embedding delay for target that is in use
	% pairRange -- range within which we've pulled source-target interactions
	% tranEntropy -- average transfer entropy
	% localTranEntropy -- local TE for each sample
	
	fprintf('%d samples in total for %d fish\n', length(S), length(unique(fileTimeAndPair(:,3))));
	
	relDistance = RelSourcePos(:,1);
	relTheta = RelSourcePos(:,2);
	% I think this is giving us the right conversions:
	if (properties.data3d)
		relPhi = RelSourcePos(:,3);
		distXY = relDistance .* cos(relPhi);
		distZ = relDistance .* sin(relPhi);
	else
		distXY = relDistance;
	end
	distInFront = distXY .* cos(relTheta); % X coordinate
	distToLeft = distXY .* sin(relTheta); % Y coordinate
	
	% Plot where all the raw positions are:
	%  Will need to turn this off when we have too many
	figure()
	% polar(relTheta, distXY, '.r'); % This is equivalent to below:
	scatter(distInFront, distToLeft, 2, localTranEntropy);
	title('Relative position of source in XY plane for target heading, coloured for TE');
	colorbar;

	% Plot the density of samples
	makePolarBinnedPlot(relTheta, distXY, ones(length(relTheta), 1), 12, 10, true, false);
	title('Density of samples in each bin (r_{XY},\theta)');
	
	% Plot the TE in XY plane
	makePolarBinnedPlot(relTheta, distXY, localTranEntropy, 12, 10, true, true);
	title('Average TE in each bin (r_{XY},\theta)');
	
	% Plot raw positions in phi-z:
	% figure()
	% polar(relPhi, relDistance, '.r');
	% title('Relative position of source in Z-phi plane for target heading');

	if (properties.data3d)
		% Plot the density of samples in distance-phi plane
		makePolarBinnedPlot(relPhi, relDistance, ones(length(relTheta), 1), 12, 10, true, false);
		title('Density of samples in each bin (r, \phi)');

		% Plot the TE in distance-phi plane
		makePolarBinnedPlot(relPhi, relDistance, localTranEntropy, 12, 10, true, true);
		title('Average TE in each bin (r, \phi)');
		xlabel('r_{XY} [mm]');
		ylabel('z [mm]');
	end
end

% Inputs:
% - thetas - angles for each sample
% - radii - radius for each sample
% - numAngleBins - how many bins to make across 2*pi
% - numRadialBins - how many bins to make up to the maximum radii
% - useMaxEntBinning - whether to make bins with approx same numbers of points (true)
%    or same size (false)
% - plotMean - if true (default) plot the mean within each bin, else plot the total (divded by area)
%    The latter is used for densities for example
function makePolarBinnedPlot(thetas, radii, valuesToPlot, numAngleBins, numRadialBins, useMaxEntBinning, plotMean)
	
	if (nargin < 6)
		useMaxEntBinning = false;
	end
	if (nargin < 7)
		plotMean = true;
	end
	
	if ((min(thetas) < -pi/2) || (max(thetas) > pi/2))
		% We're using full angular range -pi : pi
		minAngle = -pi;
		maxAngle = pi;
		extraBinForPlotWrap = true;
	else
		% We're only using -pi/2:pi/2
		minAngle = -pi/2;
		maxAngle = pi/2;
		extraBinForPlotWrap = false;
	end
	
	angleStep = (maxAngle - minAngle) / numAngleBins;
	radiusStep = max(radii) / numRadialBins;

	% Simple way to do the binning for even bins:
	% binnedAngles = floor(thetas ./ angleStep); % Gives the discrete bin for the angle
	% binnedRadii = floor(radii ./ radiusStep); % Gives the discrete bin for the radius
	
	% More general, and allowing bins to spread with points:
	if (useMaxEntBinning)
		% Space the bins for roughly same
		%  numbers of points (when examined marginally):
		sortedAngles = sort(thetas);
		binAngleEdges = [minAngle; sortedAngles(floor((1:(numAngleBins-1)).*length(sortedAngles)./numAngleBins)); maxAngle]';
		sortedRadii = sort(radii);
		binRadiusEdges = [0; sortedRadii(floor((1:(numRadialBins-1)).*length(sortedRadii)./numRadialBins)); max(radii)]';
	else
		% Space the bins equally
		binAngleEdges = minAngle:angleStep:maxAngle;
		binRadiusEdges = 0:radiusStep:max(radii);
	end
	[angleHistCounts,binnedAngles] = histc(thetas, binAngleEdges);
	[radiiHistCounts,binnedRadii] = histc(radii, binRadiusEdges);
	
	angleBinValues = 1:numAngleBins; % unique(binnedAngles);
	radiusBinValues = 1:numRadialBins; % unique(binnedRadii);
	if (extraBinForPlotWrap)
		valuesForEachBin = zeros(length(angleBinValues) + 1, length(radiusBinValues));
	else
		valuesForEachBin = zeros(length(angleBinValues), length(radiusBinValues));
	end
	numberOfSamples = 0;
	for aIndex = 1 : length(angleBinValues)
		for rIndex = 1 : length(radiusBinValues)
			indicesForThisBin = find((binnedAngles == angleBinValues(aIndex)) & (binnedRadii == radiusBinValues(rIndex)));
			if (plotMean)
				valueForThisBin = mean(valuesToPlot(indicesForThisBin));
			else
				% Plot a density: compute total then divide by area.
				valueForThisBin = sum(valuesToPlot(indicesForThisBin));
				areaOfBin = pi .* (binRadiusEdges(rIndex+1).^2 - binRadiusEdges(rIndex).^2) .* ...
						mod(abs(binAngleEdges(aIndex+1) - binAngleEdges(aIndex)), 2.*pi) ./ (2.*pi);
				valueForThisBin = valueForThisBin ./ areaOfBin;
			end
			if (length(indicesForThisBin) == 0)
				valueForThisBin = 0;
			end
			% fprintf('Mean value for r=%.1f+,theta=%.3f+ is %.3f (from %d samples)\n', binRadiusEdges(rIndex), ...
			%	binAngleEdges(aIndex), valueForThisBin, length(indicesForThisBin));
			numberOfSamples = numberOfSamples + length(indicesForThisBin);
			valuesForEachBin(aIndex, rIndex) = valueForThisBin;
		end
	end
	% Now convert these so we can plot them:
	fprintf('Found %d in total in the bins\n', numberOfSamples);
	if (extraBinForPlotWrap)
		% And add for first angle again to complete the plot
		valuesForEachBin(end,:) = valuesForEachBin(1,:);
		[THETA,RR] = meshgrid([(binAngleEdges(1:end-1)+binAngleEdges(2:end))./2, (binAngleEdges(1)+binAngleEdges(2))./2], ...
				(binRadiusEdges(1:end-1)+binRadiusEdges(2:end))./2);
	else
		[THETA,RR] = meshgrid([(binAngleEdges(1:end-1)+binAngleEdges(2:end))./2], ...
				(binRadiusEdges(1:end-1)+binRadiusEdges(2:end))./2);
	end
	[A,B] = pol2cart(THETA,RR);
	figure();
	% Old way which pinned TE values on the vertices of polygons (looks yuck)
	% surf(A,B,valuesForEachBin','edgecolor','none')
	% New way, smoothed visualisation:
	plot = pcolor(A,B,valuesForEachBin');
	plot.FaceColor = 'interp';
	set(plot, 'EdgeColor', 'none');
	xlabel('x [mm]');
	ylabel('y [mm]');
	colorbar;
	view(0,90)

end

