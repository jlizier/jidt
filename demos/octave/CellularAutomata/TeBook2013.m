% To recreate the plots in:
% chapter 5 of T. Bossomaier, L. Barnett, M. Harr&eacute; and J.T. Lizier, "An Introduction to Transfer Entropy: Information Flow in Complex Systems", Springer, to be published, 2013.
% These plots were first published in (see DirectedMeasuresChapterDemo2013.m):
% J. T. Lizier, "Measuring the dynamics of information processing on a local scale in time and space",
%  accepted for "Directed Information Measures in Neuroscience", edited by M. Wibral, R. Vicente, J. T. Lizier, to be published by Springer, 2013

clear all;

% Set up simulation options:
cells = 10000;
timeSteps = 600;
neighbourhood = 3;
caStates = 2;

% Set up options for information dynamics analysis, and which segment of the CA to plot
measureParams.k=16;
options.saveImages = true;
options.saveImagesFormat = 'pdf';
options.plotOptions.scaleColoursToSubsetOfPlot = true;
scaleColoursToExtremesDefault = false;
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault;
% Turn up the contrast so that the small values aren't disproportionately visible: (0.15, 0.30 was good, except for separable which was better with 0.15, 0.35)
options.plotOptions.scalingMainComponent = 0.15;
options.plotOptions.scalingScdryComponent = 0.30;
options.plotOptions.gammaPower = 0.5;

%%%%%%%%%
% Examining rule 54:
options.plotOptions.plotRows = 35;
options.plotOptions.plotCols = 35;
options.plotOptions.plotStartRow = 20+20;
options.plotOptions.plotStartCol = 1+10;
options.seed = 3; % Set up the random number generator to give reproducible initial states for all measurements
if (exist('initialStates/DirectedMeasuresChapterDemo2013-rule54.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/DirectedMeasuresChapterDemo2013-rule54.txt');
end
fprintf('\nStarting rule 54 ...\n');
fprintf('\nPlotting apparent transfer entropy j = 1 ...\n');
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
options.plotOptions.scalingScdryComponent = 0.35;
options.plotOptions.scalingMainComponent = 0.35; 
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
options.plotRawCa = false;
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);

% If we were going to apply this to another rule:
% fprintf('\nPress any key when ready to apply to the next rule\n')
% pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

