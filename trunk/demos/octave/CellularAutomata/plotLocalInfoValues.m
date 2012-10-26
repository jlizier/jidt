% J. Lizier, 2012.
%
% Plots the local values, both the positive and negative, in blues and reds respectively, on the current figure.
%
% Inputs:
% - localResults - local values to be plotted. Can be native octave or java array
% - plotOptions - structure (optional) containing the following variables:
%   - plotRows - how many rows to plot (default is all)
%   - plotCols - how many columns to plot (default is all)
%   - plotStartRow - which row to start plotting from (default is 1)
%   - plotStartCol - which column to start plotting from (default is 1)
%   - scaleColoursToExtremes - stretch the darkest red and blue to the max and min of the local values, 
%      regardless of how imbalanced the max and mins are. (Default is false.)
%   - mainSignVectorLength - length of the longest component (positive or negative values) for the colourmap
%   - scalingMainComponent - what proportion to make the darkest shade of the primary colour (default .25)
%   - scalingScdryComponent - what proportion to make the lighter shades with green added (default .4)
%
function plotLocalInfoValues(localResults, plotOptions)

	% Set the colormap to have blue for positive, red for negative, scaled to the max and min of our local values

	if (nargin < 2)
		plotOptions = {};
	end
	if not(isfield(plotOptions, "plotRows"))
		plotOptions.plotRows = rows(localResults);
	end
	if not(isfield(plotOptions, "plotCols"))
		plotOptions.plotCols = columns(localResults);
	end
	if not(isfield(plotOptions, "plotStartRow"))
		plotOptions.plotStartRow = 1;
	end
	if not(isfield(plotOptions, "plotStartCol"))
		plotOptions.plotStartCol = 1;
	end
	if not(isfield(plotOptions, "scaleColoursToExtremes"))
		plotOptions.scaleColoursToExtremes = false;
	end
	if not(isfield(plotOptions, "mainSignVectorLength"))
		plotOptions.mainSignVectorLength = 1024;
	end
	if not(isfield(plotOptions, "scalingMainComponent"))
		plotOptions.scalingMainComponent = .25;
	end
	if not(isfield(plotOptions, "scalingScdryComponent"))
		plotOptions.scalingScdryComponent = .4;
	end
	% Pull some options out for easier coding here:
	scalingMainComponent = plotOptions.scalingMainComponent;
	mainSignVectorLength = plotOptions.mainSignVectorLength;
	scalingScdryComponent = plotOptions.scalingScdryComponent;

	if (not(ismatrix(localResults)))
		% localResults must be a java matrix
		% Convert the local values back to octave native values, and select
		%  only the rows and cols that will be plotted (gives a big speed up over
		%  coverting all of the values)
		localResultsToPlot = javaMatrixToOctave(localResults, plotOptions.plotStartRow, plotOptions.plotStartCol, ...
				plotOptions.plotRows, plotOptions.plotCols);
		% JIDT jar file is assumed to be on the path since we already have java objects coming in
		mUtils = javaObject('infodynamics.utils.MatrixUtils');
		minLocal = mUtils.min(localResults);
		maxLocal = mUtils.max(localResults);
	else
		% Pull out the values we'll be plotting
		localResultsToPlot = localResults(plotOptions.plotStartRow:plotOptions.plotStartRow+plotOptions.plotRows - 1, ...
			plotOptions.plotStartCol:plotOptions.plotStartCol+plotOptions.plotCols - 1);
		minLocal = min(min(localResults));
		maxLocal = max(max(localResults));
	end
	printf("[max,min] for local info dynamics profile is [%.3f, %.3f]\n", maxLocal, minLocal);
	
	if (minLocal < 0)
		% We need a negative component for the colorchart
		if (plotOptions.scaleColoursToExtremes)
			% Put the darkest blues and reds at our max and min respectively
			if (abs(minLocal) > maxLocal)
				% Use a longer length for the red vector than blue
				vNegLength = mainSignVectorLength;
				vPosLength = floor(mainSignVectorLength .* maxLocal ./ abs(minLocal));
			else
				% Use a longer length for the blue vector than red
				vPosLength = mainSignVectorLength;
				vNegLength = floor(mainSignVectorLength .* abs(minLocal) ./ maxLocal);
			end
			bluePosmap = prepareColourmap(vPosLength, true, scalingMainComponent, scalingScdryComponent);
			redNegmap = flipud(prepareColourmap(vNegLength, false, scalingMainComponent, scalingScdryComponent));
			% printf("Plotting locals with scaling to extreme values (%d distinct colours for positive, %d for negative)\n", \
			%	vPosLength, vNegLength);
		else
			% Only use the darkest blue/red for which of positive or negative values
			%  had the largest absolute value. For the other, scale the colours
			%  correspondingly.
			% Initialise the full spectrum
			bluePosmap = prepareColourmap(mainSignVectorLength, true, scalingMainComponent, scalingScdryComponent);
			redNegmap = flipud(prepareColourmap(mainSignVectorLength, false, scalingMainComponent, scalingScdryComponent));
			% Then cut down the non-dominant sign's part:
			if (abs(minLocal) > maxLocal)
				% Use a longer length for the red vector than blue; chop blue down
				vPosLength = floor(mainSignVectorLength .* maxLocal ./ abs(minLocal));
				vNegLength = mainSignVectorLength;
				bluePosmap = bluePosmap(1:vPosLength,:);
			else
				% Use a longer length for the blue vector than red; chop red down
				vNegLength = floor(mainSignVectorLength .* abs(minLocal) ./ maxLocal);
				vPosLength = mainSignVectorLength;
				redNegmap = redNegmap(length(redNegmap)-vNegLength + 1:length(redNegmap),:);
			end
			% printf("Plotting locals with minor scaled to major (%d distinct colours for positive, %d for negative)\n", \
			%	vPosLength, vNegLength);
		end
		% Construct the colormap with blue for positive and red for negative
		colormap([redNegmap; bluePosmap]);
		% Now, plot the local values with the pre-prepared colormap
		imagesc(localResultsToPlot);
	else
		% We need to pin the minimum value of the blue-only plot to zero.
		bluemap = prepareColourmap(mainSignVectorLength, true, scalingMainComponent, scalingScdryComponent);
		colormap(bluemap);
		% Now, plot the local values with the pre-prepared colormap
		imagesc(localResultsToPlot, [0, maxLocal]);
	end
	colorbar
end

