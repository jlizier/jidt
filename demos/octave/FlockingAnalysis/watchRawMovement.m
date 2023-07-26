function watchRawMovement(dataFileName, properties, refreshRate, zoomIn)
% Plots a movie of raw movement data from the given input file
%
% Author: Joseph T. Lizier, 2020
%
% Inputs:
% - dataFileName - the name of the file to load
% - properties - properties object
% - refreshRate - how often to change the plotted data (default 0.1 sec)
% - zoomIn - whether to zoom into the individuals (default), or have plot take in whole field

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

if (nargin < 3)
	refreshRate = 0.1;
end

if (nargin < 4)
	zoomIn = true;
end

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
if (~zoomIn)
	fprintf('Focussing...\n');
	if (properties.data3d)
		axis([minX, maxX, minY, maxY, minZ, maxZ]);
	else
		axis([minX, maxX, minY, maxY]);
	end
	% Need to set this so that the axes don't keep updating
end
fprintf('Loading data in %s (%d missing values)\n', dataFileName, numMissing);

for t = 1:size(x, 1)
	set(gca,'NextPlot','replacechildren') ;
	if (properties.data3d)
		plot3(x(t,:), y(t,:), z(t,:), 'x')
	else
		plot(x(t,:), y(t,:), 'x')
	end
	pause(refreshRate)
end

end

