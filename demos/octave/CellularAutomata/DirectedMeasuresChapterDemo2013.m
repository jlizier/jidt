% To recreate the plots in:
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
printf('\nStarting rule 54 ...\n');
printf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
printf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
options.plotOptions.scalingScdryComponent = 0.35;
options.plotOptions.scalingMainComponent = 0.35; 
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule 18: (50, 50, 20, 900 is not too bad)
options.plotRawCa = true;
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 20;
options.plotOptions.plotStartCol = 900;
options.seed = 3; % Set up the random number generator to give reproducible initial states for all measurements
printf('\nStarting rule 18 ...\n');
printf('\nPlotting active info storage ...\n');
options.plotOptions.scalingScdryComponent = 0.45; % Make the moderately strong values easier to see:
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
printf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for complete transfer entropy j = -1 ...\n');
pause
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'transfercomplete', measureParams, options);
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value
printf('\nAll done, press any key to continue ...\n');
pause

