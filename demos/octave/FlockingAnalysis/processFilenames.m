function fileCellArray = processFilenames(fileList)
%
% Turns the fileList from the properties file (usually properties.files) into a cell array of file names. The fileList can be either:
% a. a cell array of file names, e.g.: {'file1.xlsx', 'file2.xlsx'}
% b. a call to ls or ls with an argument, e.g. ls('*.xlsx')
% c. a space or tab separated character row vector of file names
% d. a character matrix of filenames (each filename on a separate row)
%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2022, Joseph T. Lizier et al.
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

	if (iscell(fileList))
		% We're done already:
		fileCellArray = fileList;
	elseif (isvector(fileList))
		% We have a row vector of space/tab separate filenames:
		fileCellArray = strsplit(strtrim(fileList)); % extra strtrim to remove trailing \n's
	elseif (ismatrix(fileList))
		fileCellArray = {};
		for r = 1 : size(fileList, 1)
			fileCellArray{r} = strtrim(fileList(r,:));
		end
	else
		error('fileList appears to be of an incorrect format\n');
	end
end

