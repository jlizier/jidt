% function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKernel(rs)
%
% runHeartBreathRateKernel
% Version 1.0
% Joseph Lizier
% 4/8/2014
%
% Used to explore information transfer in the heart rate / breath rate example of Schreiber --
%  estimating TE using box kernel estimation (as per Schreiber).
%
% Inputs
% - rs - a scalar specifying a single, or vector specifying multiple, value of kernel widths to evaluate TE (kernel) with. Defaults to roughly match the points Schreiber used: [0.01, 0.016, 0.023, 0.032, 0.047, 0.064, 0.09, 0.12, 0.18, 0.26, 0.36, 0.50, 0.65, 1.0]
% Outputs
% - teHeartToBreath - TE (heart -> breath) for each value of r
% - teBreathToHeart - TE (breath -> heart) for each value of r


function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKernel(rs)

	tic;
	
	if (nargin < 1)
		rs = [0.01, 0.016, 0.023, 0.032, 0.047, 0.064, 0.09, 0.12, 0.18, 0.26, 0.36, 0.50, 0.65, 1.0];
	end
	
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
	
	fprintf('TE for heart rate <-> breath rate for kernel estimation with %d samples:\n', timeSteps);

	% Using a single conditional mutual information calculator is the least biased way to run this:
	teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
	teCalc.setProperty('NORMALISE', 'true'); % Normalise the individual variables. Schreiber doesn't explicitly say this is done for TE, but it is done for the raw data plots in Figure 3.
	teCalc.setProperty('DYN_CORR_EXCL', '100'); % Schreiber never mentions that dynamic correlation exclusion is used, but one suspects that it is done with 100 time steps (as per example 2), and indeed this is borne out by the results.
	% teCalc.setDebug(true);

	for rIndex = 1:length(rs)
		r = rs(rIndex);
		% Compute TE values for kernel radius r
		
		% Perform calculation for heart -> breath (lag 1)
		teCalc.initialise(1, r); % Use history length 1 (Schreiber k=1), kernel width of r normalised units
		teCalc.setObservations(heart, chestVol);
		% teHeartToBreath(rIndex) = teCalc.computeAverageLocalOfObservations();
		% Schreiber confirmed that example 2 was bias corrected (though it had no 
		%  effect for the 10000 iterations there); it seems that the same thing
		%  was done here.
		% As per the Javadocs for TransferEntropyCalculatorKernel, it is not
		%  recommended to use the bias corrected method here, since it is not clear
		%  how the corrections cancel (whereas this is done properly by the Kraskov estimator)
		teHeartToBreath(rIndex) = teCalc.computeAverageLocalOfObservationsWithCorrection();
		
		% Perform calculation for breath -> heart (lag 1)
		teCalc.initialise(1, r); % Use history length 1 (Schreiber k=1), kernel width of r normalised units
		teCalc.setObservations(chestVol, heart);
		% teBreathToHeart(rIndex) = teCalc.computeAverageLocalOfObservations();
		% Schreiber confirmed that example 2 was bias corrected (though it had no 
		%  effect for the 10000 iterations there); it seems that the same thing
		%  was done here.
		% As per the Javadocs for TransferEntropyCalculatorKernel, it is not
		%  recommended to use the bias corrected method here, since it is not clear
		%  how the corrections cancel (whereas this is done properly by the Kraskov estimator)
		teBreathToHeart(rIndex) = teCalc.computeAverageLocalOfObservationsWithCorrection();
		
		fprintf('TE(r=%.2f): heart->breath = %.3f, breath->heart = %.3f\n', r, teHeartToBreath(rIndex), teBreathToHeart(rIndex));
	end
		
	% Plot results 
	hold off;
	semilogx(rs, teHeartToBreath, '-r');
	hold on;
	% And add points we extracted from Schreiber's plot for TE(1->2):
	schreiberResultsHeartToBreath =  load('SchreiberExample3.txt');
	semilogx(schreiberResultsHeartToBreath(:,1), schreiberResultsHeartToBreath(:,2), '-b');
	semilogx(rs, teBreathToHeart, '-g');
	hold off;
	set (gca,'fontsize',26);
	xlabel('kernel width', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('TE (bits)', 'FontSize', 36, 'FontWeight', 'bold');
	axis([0.01 1 0 5]);
	h = legend('TE(heart->breath)', 'Schreiber(heart->breath)', 'TE(breath->heart)');
	set (h,'fontsize',12);
	print('heartBreathResults-kernel.eps', '-depsc');

	totaltime = toc;
	fprintf('Total runtime was %.1f sec\n', totaltime);
end

