% function plotRawCa(states, saveIt)
%
% Plot the given raw values of a cellular automata run
%
% Inputs
% - states - 2D array of states of the CA (1st index time goes along the rows, 2nd index cells go across the columns)
% - rule - Rule number (as number or a string) for the rule, used to label the file
% - plotOptions - structure (optional) containing the following variables:
%   - plotRows - how many rows to plot (default is all)
%   - plotCols - how many columns to plot (default is all)
%   - plotStartRow - which row to start plotting from (default is 1)
%   - plotStartCol - which column to start plotting from (default is 1)
% - saveIt - whether to save an eps file of the image (default false)

function plotRawCa(states, rule, plotOptions, saveIt)

	if (nargin < 2)
		plotOptions = {};
	end
	if not(isfield(plotOptions, 'plotRows'))
		plotOptions.plotRows = size(states, 1);
	end
	if not(isfield(plotOptions, 'plotCols'))
		plotOptions.plotCols = columns(states);
	end
	if not(isfield(plotOptions, 'plotStartRow'))
		plotOptions.plotStartRow = 1;
	end
	if not(isfield(plotOptions, 'plotStartCol'))
		plotOptions.plotStartCol = 1;
	end
	if (plotOptions.plotRows > 1000)
		fprintf('*** Limiting number of plotted rows to 1000\n');
		plotOptions.plotRows = 1000;
	end
	if (plotOptions.plotCols > 1000)
		fprintf('*** Limiting number of plotted columns to 1000\n');
		plotOptions.plotCols = 1000;
	end
	if (nargin < 3)
		saveIt = false;
	end

	figure(1);
	% set the colormap for black = 1, white = 0
	colormap(1 - gray(2))
	imagesc(states(plotOptions.plotStartRow:plotOptions.plotStartRow+plotOptions.plotRows - 1, ...
		plotOptions.plotStartCol:plotOptions.plotStartCol+plotOptions.plotCols - 1))
	colorbar
	fprintf('Adding colorbar to ensure that the size of the diagram matches that of local info plots\n');
	if (saveIt)
		set(gca, 'fontsize', 32);
		colorbar('fontsize', 32);
		if (ischar(rule))
			filename = sprintf('figures/raw-%s.eps', rule);
		else
			filename = sprintf('figures/raw-%d.eps', rule);
		end
		print(filename, '-depsc');
	end
end

