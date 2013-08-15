% To recreate the plots in:
% Lizier, Prokopenko and Zomaya, "A framework for the local information dynamics of distributed computation in complex systems",
%  accepted for "Guided Self-Organization: Inception", edited by M. Prokopenko, to be published by Springer, 2013

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
options.plotOptions.plotStartRow = 150;
options.plotOptions.plotStartCol = 175;
options.seed = 2;
printf('\nStarting rule 54 ...\n');
printf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'active', measureParams, options);
printf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'excess', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
options.plotOptions.scalingScdryComponent = 0.35;
options.plotOptions.scalingMainComponent = 0.35; 
measureParams.k=16; % back to default
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingMainComponent = 0.15; % Return to previous value
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'separable', measureParams, options);
printf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule 110: (60, 60, 220, 230 is good but with only one collision)
options.plotRawCa = true;
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 50+20;
options.plotOptions.plotStartCol = 800+60;
options.seed = 2;
printf('\nStarting rule 110 ...\n');
printf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'active', measureParams, options);
printf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'excess', measureParams, options);
printf('\nPress any key when ready for entropy rate ...\n');
pause
measureParams.k=16; % back to default
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'entropyrate', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'separable', measureParams, options);
printf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule phi_par:
% (75, 75, 12, 300 - is ok, but another collision would be good)
% (100, 100, 12, 500 - is very interesting, but perhaps slightly too complicated. Dynamics will lead to some strange colouring also)
options.plotRawCa = true;
options.plotOptions.plotRows = 70;
options.plotOptions.plotCols = 70;
options.plotOptions.plotStartRow = 17;
options.plotOptions.plotStartCol = 585;
options.seed = 2;
phi_par = 'feedffdec1aaeec0eef000a0e1a020a0';
phi_par_neighbourhood = 7;
phi_par_cells = 30000;
phi_par_timeSteps = 200;
measureParams.k=10; % Shorter for phi_par
printf('\nStarting rule phi_par ...\n');
printf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
printf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for complete transfer entropy j = -1 ...\n');
pause
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfercomplete', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = -3 ...\n');
pause
measureParams.j = -3;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingMainComponent = 0.35; % Make it easier to see stronger negatives
options.plotOptions.scalingScdryComponent = 0.45;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'separable', measureParams, options);
printf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value
options.plotOptions.scalingMainComponent = 0.15; % Return to previous value
measureParams.k=16; % back to default

%%%%%%%%%
% Examining rule 22:
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 150;
options.plotOptions.plotStartCol = 175;
options.seed = 2;
printf('\nStarting rule 22 ...\n');
printf('\nPlotting active info storage ...\n');
options.plotRawCa = true;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'active', measureParams, options);
printf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'excess', measureParams, options);
printf('\nPress any key when ready for entropy rate ...\n');
pause
measureParams.k=16; % back to default
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'entropyrate', measureParams, options);
printf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'transfer', measureParams, options);
printf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'separable', measureParams, options);
printf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

