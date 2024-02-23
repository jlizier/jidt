% function [data, bboxCentredData, selectedImageBoxCentred] = parseToText(filePrefix, sampleToPlot)
% 
% J Lizier - originally written 15/10/13, updated 30/9/2021.
%
% Parse the binary MNIST files, as described at http://yann.lecun.com/exdb/mnist/ into our matrix format
%
% Inputs:
% - filePrefix - (directory+)file prefix of the four label and image files (for original data is either 't10k' for test set or 'train' for training set) 
%    as downloaded from http://yann.lecun.com/exdb/mnist/ (you need to download in advance)
% - sampleToPlot - specify a particular image to plot. If not specified, a random one is chosen
%
% Outputs:
% - data - data in our format - each row contains class (1-10, as digit number + 1) in first column, followed by 28x28 pixels.
% - bboxCentredData - data in our format, as above, but with 20x20 pixels only, selected by the minimum bounding box, centred.
% - selectedImageBoxCentred - data for the sample plot

%
%  Java Information Dynamics Toolkit (JIDT)
%  Copyright (C) 2022, Joseph T. Lizier et al.
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [data, bboxCentredData, selectedImageBoxCentred] = parseToText(filePrefix, sampleToPlot)

% Constants:
bboxImageRows = 20; % Should not change this, or else we may need to worry that it goes outside the images here (which were fitted around this I think)
bboxImageColumns = 20;

% Load in the labels first:
labelsFname = [filePrefix, '-labels.idx1-ubyte'];
[labelsFid, msg] = fopen(labelsFname);
if (labelsFid == -1)
	error(['fopen error for labels file: ', msg]);
end

% Read the magic number, 2049, first: (MSB first which is big-endian)
[maginNum, count] = fread (labelsFid, 1, 'int32', 0, 'ieee-be');
if (count < 1)
	error('Magic number not read in labels file');
end
if (maginNum ~= 2049)
	error(sprintf('Magic number in labels file %d <> 2049', maginNum));
end

% Next get the number of items: (MSB first which is big-endian)
[numLabels, count] = fread (labelsFid, 1, 'int32', 0, 'ieee-be');
if (count < 1)
	error('Number of items not read in labels file');
end
fprintf('Reading %d labels in the file ...\n', numLabels);

% Now read all the labels: (unsigned byte)
[labels, count] = fread (labelsFid, numLabels, 'uint8', 0, 'ieee-be');
if (count < numLabels)
	error('Number of labels contained in labels file (%d) is smaller than reported in the file header (%d)', count, numLabels);
end

% Now open the image file
imagesFname = [filePrefix, '-images-idx3-ubyte']; % Used to be .idx3
[imagesFid, msg] = fopen(imagesFname);
if (imagesFid == -1)
	error(['fopen error for images file: ', msg]);
end

% Read the magic number, 2051, first: (MSB first which is big-endian)
[maginNum, count] = fread (imagesFid, 1, 'int32', 0, 'ieee-be');
if (count < 1)
	error('Magic number not read in images file');
end
if (maginNum ~= 2051)
	error(sprintf('Magic number in images file %d <> 2051', maginNum));
end

% Next get the number of images, rows and columns: (MSB first which is big-endian)
[headerData, count] = fread (imagesFid, 3, 'int32', 0, 'ieee-be');
if (count < 3)
	error('Number of images, rows, columns not read in labels file');
end
numImages = headerData(1);
numRows = headerData(2);
numColumns = headerData(3);
if (numImages ~= numLabels)
	error('Number of images %d said to be in images file does not match number in labels file %d', numImages, numLabels);
end
fprintf('Reading %d images of size %dx%d in the file ...\n', numImages, numRows, numColumns);

% Now read all the images: (unsigned byte).
% It seems fread fills a column before moving onto the next one, so give our desired rows and columns in reverse ...
[imageData, count] = fread (imagesFid, [numRows*numColumns, numImages], 'uint8', 0, 'ieee-be');
if (count < numImages*numRows*numColumns)
	error('Number of images contained in images file (%d) is smaller than reported in the file header (%d)', count, numImages*numRows*numColumns);
end
% Then bring our rows and columns back to our usual format:
imageData = imageData';
fprintf('Read image data with %d rows of %d columns\n', size(imageData,1), size(imageData, 2));

% Plot one image just to make sure we've got it correctly:
if (nargin < 2)
	% If we're not given a particular image to plot, pick one at random
	sampleToPlot = unidrnd(numImages);
end
fprintf('Plotting sample image %d, which is labelled %d\n', sampleToPlot, labels( sampleToPlot));
figure
imagesc(reshape(imageData(sampleToPlot,:), numRows, numColumns)');
colormap(1 - gray)
axis([0.5 numRows+0.5 0.5 numColumns+0.5]);
title(sprintf('Sample image %d',  sampleToPlot));

% Finally, reformat all the data, in our format -- classes must be indexed from 1
data = [labels + 1, imageData];
bboxCentredData = zeros(numImages, 1 + bboxImageRows*bboxImageColumns);
bboxCentredData(:,1) = data(:,1);

% Check the bounding box sizes here, and attempt bounding box centring:
xsizes = zeros(numImages, 1);
ysizes = zeros(numImages, 1);
for in = 1 : numImages
	theImage = reshape(imageData(in,:), numRows, numColumns)';
	[nonZerosX, nonZerosY] = find(theImage);
	xsizes(in) = max(nonZerosX) - min(nonZerosX) + 1; % +1 to be inclusive of the edge pixels
	ysizes(in) = max(nonZerosY) - min(nonZerosY) + 1;
	% Now do the centering in bboxImageRows x bboxImageColumns pixels:
	paddingX = bboxImageRows - xsizes(in);
	paddingY = bboxImageColumns - ysizes(in);
	% **Design decision**:
	padLessLowerIndices = false; % pad less on lower indices. Not sure whether this matters - might be different to the original data set.
	% Also, padding should not go outside theImage bounds if rows/columns set to their original value of 20.
	try
		if (padLessLowerIndices)
			x1 = min(nonZerosX) - floor(paddingX/2);
			x2 = max(nonZerosX) + ceil(paddingX/2);
			y1 = min(nonZerosY) - floor(paddingY/2);
			y2 = max(nonZerosY) + ceil(paddingY/2);
		else
			x1 = min(nonZerosX) - ceil(paddingX/2);
			x2 = max(nonZerosX) + floor(paddingX/2);
			y1 = min(nonZerosY) - ceil(paddingY/2);
			y2 = max(nonZerosY) + floor(paddingY/2);
		end
		bboxCentredImage = theImage(x1 : x2, y1 : y2);
	catch
		% Assume that we went outside the image bounds
		fprintf('Image %d bounding box failed: min(nonZerosX)=%d, max(nonZerosX)=%d, min(nonZerosY)=%d, max(nonZerosY)=%d\n', in, min(nonZerosX), max(nonZerosX), min(nonZerosY), max(nonZerosY));
		fprintf('Image %d: Pulling from x=%d:%d y=%d:%d failed - probably outside image bounds\n', in, ...
			x1, x2, y1, y2);
		fprintf('Using zeros when outside of image bounds instead.\n');
		% Use a hack to fix it:
		tempPadding = 20;
		theImage = [zeros(tempPadding, 2*tempPadding + size(theImage, 2)); ... 
			    zeros(size(theImage, 1), tempPadding), theImage, zeros(size(theImage, 1), tempPadding); ...
			    zeros(tempPadding, 2*tempPadding + size(theImage, 2))];
		bboxCentredImage = theImage(x1 + tempPadding : x2  + tempPadding, y1 + tempPadding : y2 + tempPadding);
	end
	bboxCentredData(in,2:end) = reshape(bboxCentredImage, 1, bboxImageRows*bboxImageColumns);
end
figure
hist(xsizes, 20);
title('x sizes');
figure
hist(ysizes, 20);
title('y sizes');

% And plot the bounding box centre image of the digit we plotted earlier:
figure
selectedImageBoxCentred = reshape(bboxCentredData(sampleToPlot,2:end), bboxImageRows, bboxImageColumns);
imagesc(selectedImageBoxCentred);
colormap(1 - gray)
axis([0.5 bboxImageRows+0.5 0.5 bboxImageColumns+0.5]);
title(sprintf('Sample image %d bounding box centred',  sampleToPlot));

end

