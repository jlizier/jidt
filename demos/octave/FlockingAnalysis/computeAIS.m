function ais = computeAIS(D, Dpast, properties)
% Computes active information storage from the pre-processed velocity data
%
% Author: Joseph T. Lizier, 2019
%
% Inputs:
% - D - target samples (may be multivariate as per generateObservations)
% - Dpast - target past samples (multivariate, and embedded up to k previous samples)
% - properties (required) - object with properties for the calculations,
%   with sub-members as specificied in the loadProperties.m file. If not supplied
%   the properties are loaded from loadProperties.m
%
% Outputs:
% - ais - active information storage value (MI between D and Dpast)

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

	% -- STEP 0 : create java object
	% check out http://lizier.me/joseph/software/jidt/javadocs/v1.3/
	% for the description of all classes and methods
	javaaddpath(properties.jidtJarLocation); % add JIDT path
	
	% Select correct MI class for this estimator type
	if (strcmp('kraskov', properties.estimator))
		MI_CLASS = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1';
	else
		MI_CLASS = 'infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian';
	end

	% Compute the information stored in the past of the target.
	% Creates java object of the given class 
	AIScalculator = javaObject(MI_CLASS);
	 % -- STEP 1 : set properties
	AIScalculator.setProperty('k', num2str(properties.jidt.kNNs));
	AIScalculator.setProperty('BIAS_CORRECTION', 'true'); % Used for Gaussian only
	% -- STEP 2 : initialise
	% here the parameters are the dimensionality of the series
	% in this case taken directly from the number of columns in each variable
	AIScalculator.initialise(size(D,2), size(Dpast,2));
	% -- STEP 3 : add in observations
	AIScalculator.setObservations(D, Dpast);
	% -- STEP 4 : compute the local AIS (will be bias-corrected now for either Gaussian or KSG)
	ais = AIScalculator.computeAverageLocalOfObservations(); % global (average) value
	fprintf('Mean AIS_%s (k=%d,tau=%d) = %.3f\n', properties.estimator, properties.k, properties.tau, ais);
	
	if (properties.aisNumSurrogates > 0)
		% Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations),
		%  and use this for empirical bias correction (otherwise we're relying on analytic)
		aisMeasDist = AIScalculator.computeSignificance(properties.aisNumSurrogates);
		fprintf('Null distribution: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.3f from %d surrogates)\n', ...
			aisMeasDist.getMeanOfDistribution(), aisMeasDist.getStdOfDistribution(), ...
			aisMeasDist.pValue, properties.aisNumSurrogates);
		ais = ais - aisMeasDist.getMeanOfDistribution();
		fprintf('Bias corrected Mean AIS_%s (k=%d,tau=%d) = %.3f\n', properties.estimator, properties.k, properties.tau, ais);
	end

end
