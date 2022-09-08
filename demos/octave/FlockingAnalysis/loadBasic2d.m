function [posX,posY] = loadBasic2d(dataFileName, properties)
% This script loads the raw data from a .txt file,
% preprocesses the data and save it as a .mat file.
% The txt data is assumed to have early columns with other data
%  (e.g. a timestamp in column 1), with the fish coordinates
%  starting from properties.loadBasic2d.startColumn (defaults to 2)
%  with fish 1's X and Y coordinates in the first of those columns (defaults to 2 and 3), 
%  then fish 2's X and Y coordinates in the next of those columns (defaults to 4 and 5), 
%  and so on.

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

% prompt = 'Press a key to continue';

%%% LOAD RAW DATA %%%
%fprintf('Loading data in %s\n', dataFileName); % print to check
% input(prompt);

% read the .txt as a Matlab table (assuming .txt is tab-separated)
data = load(dataFileName);
%fprintf('Size of data is %d - %d\n', ... % print table's size - %d is for int
%    size(data,1), size(data,2));         % see also %f (real) and %s (string)
% disp(data(1:5,:)); % display first 10 rows to check
% input(prompt);


%%% FORMAT DATA IN A MORE CONVENIENT WAY %%%

if (nargin == 1) || (~isfield(properties, 'loadBasic2d'))
    % Assume the first column is a datestamp
    properties.loadBasic2d.startColumn = 2;
end
startCol = properties.loadBasic2d.startColumn;
colsToSkip = startCol - 1;

numFish = (size(data,2)-colsToSkip) ./ 2; % number of fish (we know it from the raw data)
numCycles = size(data,1); % number of time steps (we also know it)
fprintf('Number of fish %d and cycles %d\n', numFish, numCycles);

% prepare variables for x and y position
% as a table [numFish x numCycles]
posX = nan(numCycles,numFish);
posY = nan(numCycles,numFish); 

% fill the position tables
for f = 1 : numFish % cycle over all fish
    % copy into new variables
    posX(:,f) = data(:,startCol+(f-1)*2);
    posY(:,f) = data(:,startCol+1+(f-1)*2);
end

% display to check
% disp(size(posX));

% disp(posX(1:5,:));
% disp(size(posY));
% disp(posY(1:5,:));
% input(prompt);

