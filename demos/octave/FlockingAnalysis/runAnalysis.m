function runAnalysis(properties)
% This high-level function generates the local transfer entropy results, first optimising
%  parameters (i.e. embedding length and delay, and source-target lag), then
%  storing local transfer entropy values for the optimised parameters.
%
% Author: Joseph T. Lizier, 2019
%
% Inputs:
% - properties (required) - object with properties for the calculations,
%   with sub-members as specificied in the loadProperties.m file. If not supplied
%   the properties are attempted to be loaded from loadProperties.m

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

if (nargin < 1)
	fprintf('No properties object supplied, attempting to load properties via a loadProperties script ...');
	% By default, just try to load properties locally
	if (exist('loadProperties') == 2)
		% there is a loadProperties script so attempt to run it to load a properties object
		loadProperties;
	else
		% there is not a loadProperties script
		error('No properties object supplied, and no loadProperties script found.');
	end
else
	% A properties argument was supplied
	if (ischar(properties))
		% We're assuming it was the name of a properties file
		if ((length(properties) > 2) && (strcmp(properties(end-1:end), '.m')))
			% Remove the '.m':
			properties(end-1:end) = [];
		end
		if (exist(properties) == 2)
			% attempt to run the properties .m file: (after making sure the properties variable is cleared; not necessary but is clean)
			propertiesFile = properties;
			clear properties;
			eval(propertiesFile);
		else
			error('%s is not an .m file we can find that can be used to load a properties object', properties);
		end
	% else
	% We assume it was the properties object.
	end
end

% Step 1: Auto-embed if required:
% Set the lag to 1 as a dummy if we are optimising over that later as well:
if (~isfield(properties, 'lag'))
	properties.lag = 1;
end
if (isfield(properties, 'kRange') || isfield(properties, 'tauRange'))
	% If any one of these two ranges weren't supplied, set the range variables
	%  to the value of the corresponding non-range variable:
	if (~isfield(properties, 'kRange'))
		properties.kRange = properties.k;
	end
	if (~isfield(properties, 'tauRange'))
		properties.tauRange = properties.tau;
	end
	% Ask generateObservations to only return the target samples for the AIS calculation
	properties.destSamplesOnly = true;
	% Optimise k and tau
	maxAIS = -inf;
	maxAISk = properties.kRange(1);
	maxAIStau = properties.tauRange(1);
	aisForKAndTau = zeros(length(properties.kRange), length(properties.tauRange));
	kIndex = 0;
	for k = properties.kRange
		kIndex = kIndex + 1;
		properties.k = k;
		minTau = min(properties.tauRange);
    	tauIndex = 0;
		for tau = properties.tauRange
			tauIndex = tauIndex + 1;
			if ((k == 1) && (tau > minTau))
				% We only need compute k=1 for a single tau
				continue;
			end
			properties.tau = tau;
			% Generate the observations for k,tau:
			[D, Dpast, ~, ~, safeDynamicCorrelationExclusionSamples] = generateObservations(properties);
			if (isempty(D))
				% There were no samples found for the given parameters, presumably k etc are too long
				continue;
			end
			% Check if we're turning on dynamic correlation exclusion:
			if (isfield(properties.jidt, 'autoDynamicCorrelationExclusion'))
				properties.jidt.dynamicCorrelationExclusion = safeDynamicCorrelationExclusionSamples;
			end
			% Compute the AIS:
			ais = computeAIS(D, Dpast, properties);
			if (ais > maxAIS)
				maxAIS = ais;
				maxAISk = k;
				maxAIStau = tau;
			end
			aisForKAndTau(kIndex, tauIndex) = ais;
		end
	end
	% Optimisation is complete:
	properties.k = maxAISk;
	properties.tau = maxAIStau;
	properties.ais = maxAIS;
	properties.destSamplesOnly = false;
	fprintf('*** Optmised k=%d and tau=%d (giving AIS=%.4f - see above for null distribution for these parameters)\n', properties.k, ...
		properties.tau, properties.ais);
else
	% Hard coded embedding parameters: compute the AIS to be saved anyway :
	% Ask generateObservations to only return the target samples for the AIS calculation
	properties.destSamplesOnly = true;
	[D, Dpast, ~, ~, safeDynamicCorrelationExclusionSamples] = generateObservations(properties);
	if (isempty(D))
		% There were no samples found for the given parameters, presumably k etc are too long
		error('No embeddings possible for the given hard coded k and tau\n');
	end
	% Check if we're turning on dynamic correlation exclusion:
	if (isfield(properties.jidt, 'autoDynamicCorrelationExclusion'))
		properties.jidt.dynamicCorrelationExclusion = safeDynamicCorrelationExclusionSamples;
	end
	% Compute the AIS:
	ais = computeAIS(D, Dpast, properties);
	properties.ais = ais;
	properties.destSamplesOnly = false;
	fprintf('*** Hard-coded values for k=%d and tau=%d to be used (giving AIS=%.4f)\n', properties.k, ...
		properties.tau, properties.ais);
end
% Check whether to continue to TE calculations:
if (~isfield(properties, 'computeAISOnly'))
	properties.computeAISOnly = false;
end
if (properties.computeAISOnly)
    fprintf('Finishing after AIS calculation only, as requested.\n');
    return;
end


% Step 2: automatically select the correct lag if required:
teNumSurrogates = properties.teNumSurrogates; % Store this for later, turn it off now
properties.teNumSurrogates = 0; % No need to run any surrogates during parameter fitting
if (isfield(properties, 'lagRange'))
	% Caller asks us to maximise the TE over a given range:
	maxTE = -inf;
	maxTElag = properties.lagRange(1);
	teForLag = zeros(length(properties.lagRange), 1);
	lagIndex = 0;
	for lag = properties.lagRange
		lagIndex = lagIndex + 1;
		properties.lag = lag;
		% Generate the observations for k,tau,lag:
		[D, Dpast, S, ~, safeDynamicCorrelationExclusionSamples] = generateObservations(properties);
		if (isempty(S))
			% There were no samples found for the given parameters, presumably k etc are too long
			continue;
		end
		% Check if we're turning on dynamic correlation exclusion:
		if (isfield(properties.jidt, 'autoDynamicCorrelationExclusion'))
			properties.jidt.dynamicCorrelationExclusion = safeDynamicCorrelationExclusionSamples;
		end
		% Compute the TE:
		te = computeTE(S, D, Dpast, properties);
		if (te > maxTE)
			maxTE = te;
			maxTElag = lag;
		end
		teForLag(lagIndex) = te;
	end
	% Optimisation is complete:
	properties.lag = maxTElag;
	properties.tranEntropy = maxTE;
	fprintf('*** Optmised lag=%d (giving TE=%.4f)\n', properties.lag, ...
		properties.tranEntropy);	
else
	fprintf('*** Hard-coded value for lag=%d to be used\n', properties.lag);
end

% 3. Compute TE with the correct parameters
% Now, once again pre-process the positional data into velocities, this time
%  saving them into the results file (by not requesting [S,D,Dpast] outputs):
generateObservations(properties);
% And load these samples (S, D, Dpast, properties, maxSourceSamplesForATarget, safeDynamicCorrelationExclusionSamples, etc) in from the saved file:
load(properties.resultsFile);
properties.teNumSurrogates = teNumSurrogates; % Allow surrogates to be computed for this final run with correct parameters
% And compute the TE again for the optimal parameters, this time
%  saving the files:
% Turn on dynamic correlation exclusion if required:
if (isfield(properties.jidt, 'autoDynamicCorrelationExclusion'))
	if (properties.jidt.autoDynamicCorrelationExclusion)
		properties.jidt.dynamicCorrelationExclusion = safeDynamicCorrelationExclusionSamples; % safeDynamicCorrelationExclusionSamples was loaded from the results file
	else
		properties.jidt.dynamicCorrelationExclusion = 0; % no dynamic correlation exclusion
	end
end
% Compute TE with no output arguments so that results are saved
computeTE(S, D, Dpast, properties);
save(properties.resultsFile, 'ais', '-append'); % Add the AIS into the results file as well
if (isfield(properties, 'kRange') || isfield(properties, 'tauRange'))
	save(properties.resultsFile, 'aisForKAndTau', '-append'); % Add the AISs computed in auto-embedding
end
if (isfield(properties, 'lagRange'))
	save(properties.resultsFile, 'teForLag', '-append'); % Add the TEs computed in optimising the source-target lag
end

end

