% function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKraskov(knns)
%
% runHeartBreathRateKraskov
% Version 1.0
% Joseph Lizier
% 22/1/2014
%
% Used to explore information transfer in the heart rate / breath rate example of Schreiber --
%  but estimates TE using Kraskov-Stoegbauer-Grassberger estimation.
%
%
% Inputs
% - kHistory - destination embedding length
% - lHistory - source embedding length
% - knns - a scalar specifying a single, or vector specifying multiple, value of K nearest neighbours to evaluate TE (Kraskov) with
% Outputs
% - teHeartToBreath - TE (heart -> breath) for each value of k nearest neighbours
% - teBreathToHeart - TE (breath -> heart) for each value of k nearest neighbours


function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKraskov(kHistory, lHistory, knns)

	starttime = tic;
	
	% Add utilities to the path
	addpath('..');

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	data = load('../../data/SFI-heartRate_breathVol_bloodOx.txt');
	
	% Restrict to the samples that Schreiber mentions:
	data = data(2350:3550,:);
	
	% Separate the data from each column:
	heart = data(:,1);
	chestVol = data(:,2);
	bloodOx = data(:,3);
	timeSteps = length(heart);
	
	fprintf('TE for heart rate <-> breath rate for Kraskov estimation with %d samples:\n', timeSteps);

	% Using a KSG estimator for TE is the least biased way to run this:
	teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');

	for knnIndex = 1:length(knns)
		knn = knns(knnIndex);
		% Compute a TE value for knn nearest neighbours
				
		% Perform calculation for heart -> breath (lag 1)
		teCalc.initialise(kHistory,1,lHistory,1,1);
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(heart), ...
					octaveToJavaDoubleArray(chestVol));
		teHeartToBreath(knnIndex) = teCalc.computeAverageLocalOfObservations();
		
		% Perform calculation for breath -> heart (lag 1)
		teCalc.initialise(kHistory,1,lHistory,1,1);
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(chestVol), ...
					octaveToJavaDoubleArray(heart));
		teBreathToHeart(knnIndex) = teCalc.computeAverageLocalOfObservations();
		
		fprintf('TE(k=%d): heart->breath = %.3f, breath->heart = %.3f nats\n', knn, teHeartToBreath(knnIndex), teBreathToHeart(knnIndex));
	end
		
	tElapsed = toc(starttime);
	fprintf('Total runtime was %.1f sec\n', tElapsed);
	
	plot(knns, teHeartToBreath, 'rx', 'markersize', 10);
	hold on;
	plot(knns, teBreathToHeart, 'mo', 'markersize', 10);
	hold off;
	legend(['TE(heart->breath)'; 'TE(breath->heart)']);
	set (gca,'fontsize',26);
	xlabel('K nearest neighbours', 'FontSize', 44, 'FontWeight', 'bold');
	ylabel('TE', 'FontSize', 44, 'FontWeight', 'bold');
	print('results.eps', '-deps', '-color');
end

