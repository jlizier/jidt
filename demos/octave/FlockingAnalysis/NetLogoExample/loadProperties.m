%
% This script loads properties for transfer entropy analysis of the data from
%  the NetLogo Flocking model.
%
% Author: Joseph T. Lizier, Sooraj Sekhar, 2019
%

%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2019, Joseph T. Lizier
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


% This script loads the default properties for the transfer entropy processing:
clear('properties');

%%%%%%%%%
% FILENAMES
%%%%%%%%%

% Input data files:
% properties.files can be:
% a. a cell array of file names, e.g.: {'file1.xlsx', 'file2.xlsx'}
% b. a call to ls or ls with an argument, e.g. ls('*.xlsx')
% c. a space or tab separated character row vector of file names
% d. a character matrix of filenames (each filename on a separate row)
properties.files = 'positions%s.txt';

% Function to read in the data files:
% loadScript must point to a function .m file that accepts two arguments
%  (the name of a file, and properties object) and returns [x,y,z] (z optional, only when 3D)
%  data where each is an array, e.g. x(time, fishIndex) indexed first by time and second by fish index.
% Use the name of the .m file after an "@" character:
properties.loadScript = @loadseparatexy;

% Is the data returned by the loadScript 3D (true) or 2D (false)?
properties.data3d = false;

% Results file - will hold the parsed velocities / relative positions, plus the
%  local transfer entropy results
properties.resultsFile = 'results.mat';

%%%%%%%%%
% PARAMETERS
%%%%%%%%%

% Distance within which to consider a pair for the info theoretic analysis (units are as per what is used in the data files)
properties.pairRange = 4; % These ones have a causal range of 3

%%%%%%%%%
% INFORMATION THEORETIC Parameters
% Only lag is used for computing lagged mutual information
% All lag, k and tau are used for transfer entropy
% k - embedding dimension of the past of the destination array.
% tau - embedding delay: time cycles separating each element in the past of the destination.
% lag - time delay between the source and target in cycles

% You can set kRange, tauRange and lagRange to ask that these are optimised by runAnalysis:
% properties.kRange = 1:10;
% properties.tauRange = 1:4; 
% properties.lagRange = 1:10;
% You can also set k, tau and lag to values that generateObservations should use
%  (although note that if this is called via runAnalysis then it will overwrite them):
properties.k = 1;
properties.tau = 1;
properties.lag = 1;

% Do we compute TE based on headings or speeds or both for the info calculations? (need at least one of these set to true)
properties.speedcalc = false;
properties.headingcalc = true;

% Additionally: Do we include the relative source position in the transfer entropy calculation (true), or
%  only the relative source heading (false)
properties.includeSourcePositionInTransfer = false;

% Do we take relative source heading and position with respect to dest heading at that same
%  time point (true, this is what we did for Crosato paper) or with respect
%  to dest heading just previous to state update (false)?
% UPDATE THIS COMMENT TO MAKE MORE SENSE. Ash and Liss say false is more biologically realistic
properties.sourceWrtSameDestTime = true;

% JIDT location:
properties.jidtJarLocation = '../../../../infodynamics.jar';

% Which estimator to use.
% Valid values are 'gaussian' (linear) or 'kraskov' (non-linear)
% properties.estimator = 'gaussian';
properties.estimator = 'kraskov';

% Properties for JIDT estimators:
properties.jidt.kNNs = 4; % Number of nearest neighbours for Kraskov algorithm: just use 4 (default)
properties.jidt.autoDynamicCorrelationExclusion = true; % Exclude nearest neighbours from at least the same target transition from being included in counts for TE. This requires the Econometrics Toolbox in Matlab in order to call autocorr. Set to false if you don't have that to avoid this being called.

properties.aisNumSurrogates = 0; % Number of surrogate calculations to run for AIS (just to see the noise floor. 0 means skip)
properties.teNumSurrogates = 0; % Number of surrogate calculations to run for TE (just to see the noise floor. 0 means skip)

% Exit the analysis after calculating the AIS (i.e. don't compute TE)
properties.computeAISOnly = false;
