function watchLocalTEs(dataFileIndex, properties, refreshRate, zoomIn, plotTEIn)
% Plots a movie of local TEs in the swarm from the given input file
%
% Author: Joseph T. Lizier, 2020
%
% Inputs:
% - dataFileIndex - the index of the file to load, from the list of files listed in the properties file. Default 1.
% - properties - object with properties for the calculations,
%   with sub-members as specificied in the loadProperties.m file. If not supplied
%   the properties are loaded from loadProperties.m
% - refreshRate - how often to change the plotted data (default 0.1 sec)
% - zoomIn - whether to zoom into the individuals (default), or have plot take in whole field
% - plotTEIn - whether to plot the average TE into a target (true, default) or average TE out from a source (false)

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
		dataFileIndex = 1;
	end

	if (nargin < 2)
		fprintf('No properties object supplied, attempting to load properties via a loadProperties script ...\n');
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

	if (nargin < 3)
		refreshRate = 0.1;
	end

	% Call utility to put filename lists in a common format
	files = processFilenames(properties.files);
	dataFileName = files{dataFileIndex};

	% load preprocessed data using the function specified in properties.loadScript
	if (properties.data3d)
		[x,y,z] = feval(properties.loadScript, dataFileName, properties);
		numMissing = sum(sum(isnan([x,y,z])));
		maxZ = max(z(:));
		minZ = min(z(:));
	else
		[x,y] = feval(properties.loadScript, dataFileName, properties);
		numMissing = sum(sum(isnan([x,y])));
	end
	maxX = max(x(:));
	minX = min(x(:));
	maxY = max(y(:));
	minY = min(y(:));
	if (nargin < 4)
		zoomIn = true;
	end
	figure()
	if (~zoomIn)
		fprintf('Focussing...\n');
		if (properties.data3d)
			axis([minX, maxX, minY, maxY, minZ, maxZ]);
		else
			axis([minX, maxX, minY, maxY]);
		end
		% Need to set this so that the axes don't keep updating
	end
	fprintf('%d missing values\n', numMissing);

	if (nargin < 5)
		plotTEIn = true;
	end

	% Work out the range of TEs for this data file:
	teMin = min(localTranEntropy(find(fileTimeAndPair(:,1) == dataFileIndex)));
	teMax = max(localTranEntropy(find(fileTimeAndPair(:,1) == dataFileIndex)));

	% Loop over all the time steps in this data file
	numFish = size(x,2);
	cb = colorbar;
    xlabel('x');
    ylabel('y');
    cb.Label.String = 'Av local TE';
	% caxis([teMin teMax]); % These are likely too extreme for the averages
	for t = 1:size(x, 1)
		set(gca,'NextPlot','replacechildren') ;
		% Now loop over all fish as either source or target:
		averageTEs = zeros(1,numFish);
		for f = 1:numFish
			% Find which interactions involve TE into or out from this fish at this time step
			if (plotTEIn)
				% For TE in, match the target:
				rowIDs = find((fileTimeAndPair(:,1) == dataFileIndex) & (fileTimeAndPair(:,2) == t) & ...
						(fileTimeAndPair(:,3) == f));
			else
				% For TE out, match the source:
				rowIDs = find((fileTimeAndPair(:,1) == dataFileIndex) & (fileTimeAndPair(:,2) == t) & ...
						(fileTimeAndPair(:,4) == f));
			end
			% Now average the TE into or out of this fish:
			averageTEs(f) = mean(localTranEntropy(rowIDs));
		end
		if (properties.data3d)
			scatter3(x(t,:), y(t,:), z(t,:), 5, averageTEs);
		else
			scatter(x(t,:), y(t,:), 5, averageTEs)
		end
        if (plotTEIn)
            title(sprintf('Average TEs into each individual at time %d, coloured for TE', t));
        else
            title(sprintf('Average TEs out from each individual at time %d, coloured for TE', t));
        end
		pause(refreshRate)
	end

end

