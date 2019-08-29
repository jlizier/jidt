function [x,y] = loadseparatexy(filename, properties)
% loadseparatexy loads 2D fish data from 2 separate txt files (one for x, one for y)
%  where in each file time increases down the
%  rows and then across the columns we have position columns for each 
%  fish in turn, i.e. in position x file:
%  <fish1x>, <fish2x>, <fish3x>, etc
% 
% Inputs:
%  - filename - the name template of the file to load, with %s where 'x' and 'y' should be filled in
%  - properties (not required) - properties object (may be required for other file loaders)
% Outputs:
%  - x - 2D array, each row contains x position for each fish (in columns)
%  - y - as per x

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

	xfilename = sprintf(filename, 'x');
	x = load(xfilename);

	yfilename = sprintf(filename, 'y');
	y = load(yfilename);

end

