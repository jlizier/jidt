% Launch the AutoAnalyser from within Matlab. You will need to have Matlab open at the demos/AutoAnalyser folder.
% This is useful where you don't have a separate Java Runtime Environment installed, and so utilise Matlab's

warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('../../infodynamics.jar');
autoAnalyser = javaObject('infodynamics.demos.autoanalysis.AutoAnalyserLauncher', false);
fprintf('Be warned - closing the AutoAnalyser applet seems to close Matlab completely, so only do so when ready to exit Matlab!\n');

