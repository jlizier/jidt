% This script downloads all of the MNIST data and parses it into a nice
% form for us to read, with all training and test data in one file

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

fprintf('Downloading original MNIST data sets...\n');

downloadPrefix = 'http://yann.lecun.com/exdb/mnist/';
files = {'train-images-idx3-ubyte.gz', 'train-labels-idx1-ubyte.gz', ...
    't10k-images-idx3-ubyte.gz', 't10k-labels-idx1-ubyte.gz'};

for ix = 1:length(files)
    out_fname = websave(files{ix},[downloadPrefix, files{ix}]);
end

fprintf('Downloads complete, now unzipping ...\n');

for ix = 1:length(files)
    gunzip(files{ix});

end

fprintf('Unzipping complete, now parsing ...\n');

[trainingData, trainingBboxCentredData, selectedImageBoxCentred] = parseToMatrix('train', 1);
[testData, testBboxCentredData, selectedImageBoxCentred] = parseToMatrix('t10k', 1);

fprintf('Parsing complete, now saving ...\n');

trialAndTestData = [trainingData; testData];
save('trialAndTest-processedData.mat', 'trialAndTestData');
trialAndTestData = [trainingBboxCentredData; testBboxCentredData];
save('trialAndTest-processedData-bboxCentred.mat', 'trialAndTestData');

fprintf('\nSaving complete, all done.\n');
fprintf('See trialAndTest-processedData.mat and trialAndTest-processedData-bboxCentred.mat for the final files\n');
fprintf('You can delete the original files ending in .gz and -ubyte if you wish\n');

