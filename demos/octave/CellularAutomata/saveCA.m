%
% Utility function to save raw CA data to a text file format.
% Unfortunately MAtlab's -ascii format only writes full doubles, making much larger files than we need.
% So I wrote this.

function saveCA(filename, caStates)

	warning('off','MATLAB:Java:DuplicateClass');
	javaaddpath('../../../infodynamics.jar');
	addpath('..');

	arrayWriter = javaObject('infodynamics.utils.ArrayFileWriter');
	arrayWriter.makeIntMatrixFile(octaveToJavaIntMatrix(caStates), filename);
	
end

