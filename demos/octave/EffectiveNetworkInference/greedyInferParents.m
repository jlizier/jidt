%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Infer the parent source variables to a given target, using the greedy/iterative/multivariate algorithm with TE.
% This a simplistic implementation of the full algorithm implemented in IDTxl - https://github.com/pwollstadt/idtxl -
%  you are referred to IDTxl for an implementation with full features available.
%
% Inputs:
% - data: multivariate data, indexed by time,variableNumber, and possibly trialNumber.
% - parameters: an object containing the expected properties, or a string
%    describing the filename to run load this object in. Can include:
%      - parameters.calcType: which estimator type to use, from options
%      'discrete', 'gaussian', 'ksg', 'granger' --
%      granger is equivalent to gaussian, but performed with Oliver Cliff's
%      toolkit, to include proper autocorrelation correction.
%      - parameters.timePointsToSkipAtStart: number of time points we'll skip
%      at the start (default 0)
%      - parameters.timePointsToSkipAtEnd: number of time points we'll skip
%      at the end (default 0)
%      - parameters.numDiscreteBins: alphabet size for the discrete variables when used. (default 2)
%      - parameters.k: target embedding length, can be 'auto' to indicate auto-embedding (default 1)
%      - parameters.k_max: max target embedding length to use when parameters.k == 'auto' (default 10)
%      - parameters.numSurrogates: number of surrogates to run, or 0 for analytic surrogates (default 1000)
%      - parameters.maxDynCorrExclLags: maximum length of dynamic correlation exclusion, which will be auto-fitted (default 50)
%      - parameters.jidtLocation: path to the JIDT folder 
%      - parameters.gcToolkitLocation: Oliver's toolkit location for Granger
%      - parameters.debug: whether to print debugging results as parents are inferred (default true)
% - targetIndex: which column number to run the inference for
% - uncorrectedThresholdForOneTarget: p-value threshold (where 0 is most significant) to select sources. We will Bonferroni correct this over sources here (but user should correct over targets if they wish)

function [parentSet, results, pValues, otherResults] = greedyInferParents(data, parameters, targetIndex, uncorrectedThresholdForOneTarget)

tic

if ischar(parameters)
    % Assume that this string contains a filename which when run will load
    % a properties object for this run
    eval(['run ', parameters]);
end

% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
%javaaddpath('/home/joseph/temp/jidt-master/jidt/infodynamics.jar');
javaaddpath([parameters.jidtLocation, 'infodynamics.jar']);
% Add utilities to the path
% addpath('/home/joseph/temp/jidt-master/jidt/demos/octave');
addpath([parameters.jidtLocation, 'demos/octave']);

% Add Oliver Cliff's toolkit to path
addpath(genpath(parameters.gcToolkitLocation));

% Set parameter defaults:
if ~isfield(parameters, 'numDiscreteBins')
    parameters.numDiscreteBins = 2;
end
if ~isfield(parameters, 'k')
    parameters.k = 1;
end
if ~isfield(parameters, 'k_max')
    parameters.k_max = 10;
end
if ~isfield(parameters, 'numSurrogates')
    parameters.numSurrogates = 1000;
end
if ~isfield(parameters, 'verbose')
    parameters.verbose = true;
end
if ~isfield(parameters, 'timePointsToSkipAtStart')
    parameters.timePointsToSkipAtStart = 0;
end
if ~isfield(parameters, 'timePointsToSkipAtEnd')
    parameters.timePointsToSkipAtEnd = 0;
end
if ~isfield(parameters, 'maxDynCorrExclLags')
    parameters.maxDynCorrExclLags = 50;
end

T = size(data,1); % TIMEPOINTS
N = size(data,2); % NODES

hasMultipleTrials = (length(size(data)) > 2);
if (hasMultipleTrials)
    R = size(data,3); % TRIALS
    % Analyse data(time, variables), but skip the first
    % and last few steps
    data = data(parameters.timePointsToSkipAtStart+1:end-parameters.timePointsToSkipAtEnd,:,:);
else
    % Analyse data(time, variables), but skip the first
    % and last few steps
    data = data(parameters.timePointsToSkipAtStart+1:end-parameters.timePointsToSkipAtEnd,:);
end

threshold = uncorrectedThresholdForOneTarget / (N-1); % Bonferroni correcting the threshold
parentSet = [];

pValues = [];
results = [];
otherResults.k_history = 1;

candidates = [1:targetIndex-1, targetIndex+1:N]; % All sources except the targets are candidates as parents

% Set boolean flags for which calculator we are doing
is_jidt = true;
is_discrete = false;
is_ksg= false;
javaConverterSingleArray = 'octaveToJavaDoubleArray';
javaConverterMatrix = 'octaveToJavaDoubleMatrix';
if (strcmp(parameters.calcType, 'granger'))
    is_jidt = false;
    if (~ischar(parameters.k))
    % if (ischar(parameters.k) && ~strcmp(parameters.k, 'auto')) % I don't think this logic was correct
        parameters.k = char(string(parameters.k));
    end
    if (hasMultipleTrials)
        error('Granger calculator does not currently support multiple trials');
    end
elseif (strcmp(parameters.calcType, 'discrete'))
    is_discrete = true;
    javaConverterSingleArray = 'octaveToJavaIntArray';
    javaConverterMatrix = 'octaveToJavaIntMatrix';
elseif (strcmp(parameters.calcType, 'ksg'))
    is_ksg = true;
end

acfDecayTimes = -1 * ones(N,1);
% Grab the ACF decay time for the target
if (is_jidt)
    if (~hasMultipleTrials)
        acfDecayTimes(targetIndex) = computeAcfDecayTime(data(:, targetIndex), parameters);
    else
        acfDecayTimes(targetIndex) = computeAcfDecayTime(squeeze(data(:, targetIndex, :)), parameters);
    end
    k_history = parameters.k;
end

if (parameters.verbose)
    fprintf('Beginning greedy selection of parents for %d with threshold %.6f\n', targetIndex, threshold);
end

% LOOP 1 -- iterating over rounds of source selection
while ~isempty(candidates)
    % Whilst there are other candidates left (and we haven't quit)
    
    % Construct the calculator and set properties:
    if (is_jidt)
        if (isempty(parentSet))
                % Just doing pairwise TEs this round
                if (is_discrete)
                    if (ischar(k_history)) % assume is 'auto'
                        error('Autoembedding not supported for discrete calculator at the moment');
                    end
                    calc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', parameters.numDiscreteBins, k_history);
                else
                    if (is_ksg)
                        calc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
                    else
                        calc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian');
                    end
                    if (ischar(k_history)) % assume is 'auto'
                        calc.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS_DEST_ONLY');
                        calc.setProperty('AUTO_EMBED_K_SEARCH_MAX', string(parameters.k_max));
                    else
                        calc.setProperty('k_HISTORY', string(k_history));
                    end                        
                end
        else
                % We're doing conditional TEs this round, conditioned on
                % length(parentSet) other sources
                if (is_discrete)
                    calc = javaObject('infodynamics.measures.discrete.ConditionalTransferEntropyCalculatorDiscrete', parameters.numDiscreteBins, k_history, length(parentSet));
                else
                    if (is_ksg)
                        calc = javaObject('infodynamics.measures.continuous.kraskov.ConditionalTransferEntropyCalculatorKraskov');
                    else
                        calc = javaObject('infodynamics.measures.continuous.gaussian.ConditionalTransferEntropyCalculatorGaussian');
                    end
                    % Assume we have saved the relevant k after the first pairwise calculation
                    calc.setProperty('k_HISTORY', string(k_history));
                    % Set up the correct number of conditionals
                    calc.setProperty(calc.COND_EMBED_LENGTHS_PROP_NAME, strjoin(string(ones(length(parentSet), 1)), ','));
                    calc.setProperty(calc.COND_EMBED_DELAYS_PROP_NAME, strjoin(string(ones(length(parentSet), 1)), ','));
                    calc.setProperty(calc.COND_DELAYS_PROP_NAME, strjoin(string(ones(length(parentSet), 1)), ','));
                end
        end
        if (~hasMultipleTrials)
            destination = feval(javaConverterSingleArray, data(:, targetIndex));
            conditionals = feval(javaConverterMatrix, data(:, parentSet));
        end
    else
        destination = data(:, targetIndex);
        conditionals = data(:, parentSet);
    end
        
    thisRoundTEResults = zeros(1, length(candidates));
    thisRoundpValResults = zeros(1, length(candidates)); % Only used for Granger
    
    % LOOP 2 -- checking (conditional) TE from all current candidates,
    % given current parent set
    for sIndexInCandidates = 1:length(candidates)
        % For each candidate:

        if is_jidt
            if (acfDecayTimes(candidates(sIndexInCandidates)) < 0)
                % We haven't computed the ACF decay time for this source yet
                if (~hasMultipleTrials)
                    acfDecayTimes(candidates(sIndexInCandidates)) = ...
                        computeAcfDecayTime(data(:,candidates(sIndexInCandidates)), parameters);
                else
                    acfDecayTimes(candidates(sIndexInCandidates)) = ...
                        computeAcfDecayTime(squeeze(data(:,candidates(sIndexInCandidates),:)), parameters);
                end
            end
            if is_ksg
                calc.setProperty('DYN_CORR_EXCL', num2str(max(acfDecayTimes([targetIndex,candidates(sIndexInCandidates),parentSet]))));
            end
            % 3. Initialise the calculator for (re-)use:
            calc.initialise();
            % 4. Supply the sample data:
            calc.setDebug(true);
            if (~hasMultipleTrials)
                sourceTimeSeries = feval(javaConverterSingleArray, data(:, candidates(sIndexInCandidates)));
                if (isempty(parentSet))
                    calc.setObservations(sourceTimeSeries, destination);
                else
                    calc.setObservations(sourceTimeSeries, destination, conditionals);
                end
            else
                calc.startAddObservations();
                for numTrial = 1 : R
                    destination = feval(javaConverterSingleArray, squeeze(data(:, targetIndex, numTrial)));
                    conditionals = feval(javaConverterMatrix, squeeze(data(:, parentSet, numTrial)));
                    sourceTimeSeries = feval(javaConverterSingleArray, squeeze(data(:, candidates(sIndexInCandidates),numTrial)));
                    if (isempty(parentSet))
                        calc.addObservations(sourceTimeSeries, destination);
                    else
                        calc.addObservations(sourceTimeSeries, destination, conditionals);
                    end
                end
                calc.finaliseAddObservations();
            end
            calc.setDebug(false);
            % 5. Compute the estimate:
            result = calc.computeAverageLocalOfObservations();
            thisRoundTEResults(sIndexInCandidates) = result;
            
            if ischar(k_history)
                k_history = calc.getProperty('k_HISTORY');
                % We autoembedded the target history if we were going to - now grab the
                % determined value to use next time
                calc.setProperty('AUTO_EMBED_METHOD', 'NONE');
                fprintf('Target history embedding set to %s\n', k_history);
            end
            otherResults.k_history = k_history;

        else
            % Compute Granger via Oliver's toolkit:
            sourceTimeSeries = data(:, candidates(sIndexInCandidates));
            [result,pval] = mvgc(destination,sourceTimeSeries,conditionals, ...
                'p',parameters.k,'q','1','test','modified','surrogates',parameters.numSurrogates);
            thisRoundTEResults(sIndexInCandidates) = result;
            thisRoundpValResults(sIndexInCandidates) = pval;
            % TODO need to readout the k history here
        end        
    end
    
    % Find the strongest source out of these candidates:
    if is_jidt
        [maxTE, maxIndex] = max(thisRoundTEResults);
    else
        % pval is higher for more significant. We will use this to
        % determine the source selection, since it corrects the raw
        % measure values for autocorrelation in this toolkit
        [pValue, maxIndex] = max(thisRoundpValResults);
        maxTE = thisRoundTEResults(maxIndex);
    end
    strongestSource = candidates(maxIndex);
    % fprintf('Strongest source is %d with TE %.5f\n', strongestSource, maxTE);
    % Check significance of this source (first need to set up its calculator again):
    if is_jidt
        calc.initialise();
        if (~hasMultipleTrials)
            sourceTimeSeries = feval(javaConverterSingleArray, data(:, strongestSource));
            if (isempty(parentSet))
                calc.setObservations(sourceTimeSeries, destination);
            else
                calc.setObservations(sourceTimeSeries, destination, conditionals);
            end
        else
            calc.startAddObservations();
            for numTrial = 1 : R
                destination = feval(javaConverterSingleArray, squeeze(data(:, targetIndex, numTrial)));
                conditionals = feval(javaConverterMatrix, squeeze(data(:, parentSet, numTrial)));
                sourceTimeSeries = feval(javaConverterSingleArray, squeeze(data(:, strongestSource, numTrial)));
                if (isempty(parentSet))
                    calc.addObservations(sourceTimeSeries, destination);
                else
                    calc.addObservations(sourceTimeSeries, destination, conditionals);
                end
            end
            calc.finaliseAddObservations();
        end
        maxTE = calc.computeAverageLocalOfObservations();
        if (parameters.numSurrogates == 0)
            measDist = calc.computeSignificance();
        else
            measDist = calc.computeSignificance(parameters.numSurrogates);
        end
        pValue = 1 - measDist.pValue; % Complementing the p value so it's the proportion of null the measure is greater than
    end
    
    if (pValue > 1 - threshold)
        % We add this source to the parent set
        if (parameters.verbose)
            fprintf('Selected source %d, with TE(%d->%d | %s)=%.5f, p-value=%.5f (conditioning on %d parents)\n', ...
                strongestSource, strongestSource, targetIndex, strjoin(string(parentSet)), maxTE, pValue, length(parentSet));
        end
        candidates(maxIndex) = []; % Remove this source from the candidates
        parentSet = [parentSet, strongestSource];
        results = [results, maxTE];
        pValues = [pValues, pValue];
    else
        % Source was not significant, so we quit
        if (parameters.verbose)
            fprintf('-- Max TE source %d was not significant (TE(%d->%d | %s)=%.5f, p-value=%.6f (threshold %.6f)), quitting\n', ...
                strongestSource, strongestSource, targetIndex, strjoin(string(parentSet)), maxTE, pValue, 1-threshold);
        end
        break;
    end
    
    toc
end

if (parameters.verbose)
    fprintf('\nFinal selected parents: %s -> %d\n', strjoin(string(parentSet)), targetIndex);
end

end

% Returns the first time the ACF dips below 1/e for the series x,
% or if there are multiple series for x we take the mean across all of them
function acfDecayTime = computeAcfDecayTime(x, parameters)
    numTrials = size(x,2);
    acfDecayTimes = zeros(1,numTrials);
    
    for trial = 1 : numTrials
        [acfValues, ~] = autocorr(x(:, trial), 'NumLags', parameters.maxDynCorrExclLags);
        acfDecayTimes(trial) = parameters.maxDynCorrExclLags; % Default is max value
        for t = 1 : parameters.maxDynCorrExclLags
            if (acfValues(t) < exp(-1))
                acfDecayTimes(trial) = t;
                break;
            end
        end
    end
    acfDecayTime = round(mean(acfDecayTimes));
end
