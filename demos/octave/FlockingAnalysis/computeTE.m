function tranEntropy = computeTE(S, D, Dpast, properties)
% Computes transfer entropy from the pre-processed velocity data
%
% Author: Emanuele Crosato, Joseph T. Lizier, 2019
%
% Inputs:
% - S - source samples (may be multivariate as per generate3DObservations)
% - D - target samples (may be multivariate as per generate3DObservations)
% - Dpast - target past samples (multivariate, and embedded up to k previous samples)
% - properties (required) - object with properties for the calculations,
%   with sub-members as specificied in the loadProperties.m file. If not supplied
%   the properties are loaded from loadProperties.m
%
% Outputs:
% - te - transfer entropy value (conditional MI from S (maybe plus RelSourcePos) to D given Dpast). If not requested, then
%   the te (and an array of local values) is saved to properties.resultsFile

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
	
	% the following is the java class for conditional mutual information in JIDT
	% transfer entropy is infact mutual information conditioning on the past
	% of the destination
	if (strcmp('kraskov', properties.estimator))
		CMI_CLASS = 'infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov1';
	else
		CMI_CLASS = 'infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian';
	end
	
	TEcalculator = javaObject(CMI_CLASS); % creates a java object of the given class 
	 % -- STEP 1 : set properties
	TEcalculator.setProperty('k', num2str(properties.jidt.kNNs));
	TEcalculator.setProperty('BIAS_CORRECTION', 'true'); % Used for Gaussian only
	if (isfield(properties.jidt, 'dynamicCorrelationExclusion'))
		% We'll ensure samples from the same target transition aren't included in nearest neighbour counts
		%  (it will exclude some others as well, but this only adds some small noise to the calculation)
		TEcalculator.setProperty('DYN_CORR_EXCL', num2str(properties.jidt.dynamicCorrelationExclusion));
	end
	% -- STEP 2 : initialise
	% here the parameters are the dimensionality of the series
	% in this case taken directly from the number of columns in each variable
	TEcalculator.initialise(size(S,2), size(D,2), size(Dpast,2));
	% -- STEP 3 : add in observations
	TEcalculator.setObservations(S, D, Dpast);
	% -- STEP 4 : compute the local entropies
	tranEntropy = TEcalculator.computeAverageLocalOfObservations(); % global (average) value
	
	fprintf('Mean TE_%s (k=%d,tau=%d,lag=%d) = %.4f\n', ...
		properties.estimator, properties.k, properties.tau, properties.lag, tranEntropy);

	if (properties.teNumSurrogates > 0)
		% Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
		measDist = TEcalculator.computeSignificance(properties.teNumSurrogates);
		fprintf('Null distribution: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)\n', ...
			measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), ...
			measDist.pValue, properties.teNumSurrogates);
		pValue = measDist.pValue;
		meanOfSurrogates = measDist.getMeanOfDistribution();
		stdOfSurrogates = measDist.getStdOfDistribution();
	else
		pValue = 1;
		meanOfSurrogates = 0;
		stdOfSurrogates = 0;
	end

	if (nargout >= 1)
		% Supply the samples back to the caller
		%  (the caller is probably trying to optimise parameters at the moment)
		% Nothing to do then actually...
	else
		% We're going to save the results instead
		
		% First generate the local values to save as well:
		localTranEntropy = TEcalculator.computeLocalOfPreviousObservations(); % local values
		
		% save results
		save(properties.resultsFile, 'tranEntropy', 'localTranEntropy', 'pValue', 'meanOfSurrogates', 'stdOfSurrogates', '-append');
		fprintf('Transfer entropy saved in %s\n', properties.resultsFile);
	end

end
