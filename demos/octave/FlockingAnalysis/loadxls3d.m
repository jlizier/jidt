function [x,y,z] = loadxls3d(filename, properties)
% loadxls3d loads 3D fish data from an xls file where time increases down the
%  rows and then across the columns we have 3 x,y,z position columns for each 
%  fish in turn, i.e.:
%  <timestamp1>, <fish1x>, <fish1y>, <fish1z>, <fish2x>, <fish2y>, <fish2z>, etc
%  <timestamp2>, <fish1x>, <fish1y>, <fish1z>, <fish2x>, <fish2y>, <fish2z>, etc
% 
% Inputs:
%  - filename - the name of the file to load
%  - properties (not required) - properties object (may be required for other file loaders)
% Outputs:
%  - x - 2D array, each row contains x position for each fish (in columns)
%  - y - as per x
%  - z - as per x

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

	[data,txt,raw] = xlsread(filename);

	% Make sure we preprocess data to have x,y, and z as 2D arrays of x(timeStep, fishID):
	
	% timeSteps = size(data,1); % number of rows; not required
	xFishIndex = 2 : 3 : size(data,2); % which indices are the x values for different fish
	x = data(:,xFishIndex);
	y = data(:,xFishIndex+1);
	z = data(:,xFishIndex+2);

end

