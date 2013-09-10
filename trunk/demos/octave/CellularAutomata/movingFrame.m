% To recreate the plots in:
% Lizier and Mahoney, "Moving Frames of Reference, Relativity and Invariance in Transfer Entropy and Information Dynamics",
%  submitted to Entropy, 2012

clear all;

% Set up options for TE j=-1 cell to the right, and which segment of the CA to plot
measureParams.k=16;
measureParams.j = -1;
options.plotOptions.plotRows = 60;
options.plotOptions.plotCols = 60;
options.plotOptions.plotStartRow = 125;
options.plotOptions.plotStartCol = 125;
options.seed = 1;
if (exist('initialStates/MovingFrameDemo2013-initialState.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/MovingFrameDemo2013-initialState.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
options.saveImages = true;
options.plotOptions.scaleColoursToSubsetOfPlot = true;
% Turn up the contrast so that the small values aren't disproportionately visible:
options.plotOptions.scalingMainComponent = 0.1;
options.plotOptions.scalingScdryComponent = 0.15;
options.movingFrameSpeed = 0;

% Examining rule 54, stationary frame:
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'active', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j=1 channel...\n');
pause
measureParams.j = 1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
fprintf('\nCopy the figures from the figures directory, then press any key')
fprintf('\nwhen ready for the measures in a moving frame of reference\n')
pause

options.movingFrameSpeed = 1;
% Examining rule 54, f=1 frame:
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'active', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
% Also plot TE for channel j=0 in the moving frame:
measureParams.j = 0;
fprintf('\nPress any key when ready for apparent transfer entropy with j=0...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy with j=0...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);

