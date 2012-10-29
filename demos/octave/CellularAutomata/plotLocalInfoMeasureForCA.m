% function plotLocalInfoMeasureForCA(neighbourhood, base, rule, cells, timeSteps, measureId, measureParams, options)
% 
% Plot one run of the given CA and compute and plot a local information dynamics profile for it
%
% Inputs:
% - neighbourhood - neighbourhood size for the rule (ECA has neighbourhood 3).
%     For an even size neighbourhood (meaning a different number of neighbours on each side of the cell),
%     we take an extra cell from the lower cell indices (i.e. from the left).
% - base - number of discrete states for each cell (for binary states this is 2)
% - rule - supplied as either:
%     a. an integer rule number if <= 2^31 - 1 (Wolfram style; e.g. 110, 54 are the complex ECA rules)
%     b. a HEX string, e.g. phi_par from Mitchell et al. is "0xfeedffdec1aaeec0eef000a0e1a020a0" (note: the leading 0x is not required)
% - cells - number of cells in the CA
% - timeSteps - number of rows to execute the CA for (including the random initial row)
% - measureId - which local info dynamics measure to plot - can be a string or an integer as follows:
%    - "active", 0 - active information storage (requires measureParams.k)
%    - "transfer", 1 - apparent transfer entropy (requires measureParams.k and j)
%    - "transfercomplete", 2 - complete transfer entropy (requires measureParams.k and j)
%    - "separable", 3 - separable information (requires measureParams.k)
%    - "all", -1 - plot all measures
% - measureParams - a structure containing options as described for each measure above:
%    - measureParams.k - history length for information dynamics measures
%    - measureParams.j - we measure information transfer across j cells to the right per time step
% - options - a stucture containing a range of other options, i.e.:
%    - plotOptions - structure as defined for the plotRawCa function
%    - seed - state for the random number generator used to set the initial condition of the CA (use this
%       for reproducibility of plots, or to produce profiles for several different measures of the same CA raw states).
%       We set rand("state", options.seed) if options.seed is supplied, and restore the previous seed afterwards.
%    - plotRawCa - default true
%    - saveImages - whether to save the plots or not (default false)
%    - movingFrameSpeed - to investigate a moving frame of reference (default 0) (as in Lizier & Mahoney paper)

function plotLocalInfoMeasureForCA(neighbourhood, base, rule, cells, timeSteps, measureId, measureParams, options)

	tic

	if (nargin < 8)
		options = {};
	end
	if not(isfield(options, "plotOptions"))
		options.plotOptions = {}; % Create it ready for plotRawCa etc
	end
	if not(isfield(options, "saveImages"))
		options.saveImages = false;
	end
	if not(isfield(options, "plotRawCa"))
		options.plotRawCa = true;
	end

	% Add utilities to the path
	addpath("..");

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	% Simulate and plot the CA
	if (isfield(options, "seed"))
		previousSeed = rand("state");
		caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false, options.seed);
		rand("state", previousSeed);
	else
		caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false);
	end
	if (options.plotRawCa)
		plotRawCa(caStates, options.plotOptions, options.saveImages);
	end
	figNum = 2;
	toc

	% convert the states to a format usable by java:
	caStatesJInts = octaveToJavaIntMatrix(caStates);
	toc
	
	plottedOne = false;

	% Make the local information dynamics measurement
	if ((ischar(measureId) && (strcmpi("active", measureId) || strcmpi("all", measureId))) || \
	    ((measureId == 0) || (measureId == -1)))
		% Compute active information storage
		activeCalc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculator', base, measureParams.k);
		activeCalc.initialise();
		activeCalc.addObservations(caStatesJInts);
		avActive = activeCalc.computeAverageLocalOfObservations();
		printf("Average active information storage = %.4f\n", avActive);
		javaLocalValues = activeCalc.computeLocalFromPreviousObservations(caStatesJInts);
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(javaLocalValues, options.plotOptions);
		if (options.saveImages)
			print("figures/active.eps", "-color", "-deps");
		end
		plottedOne = true;
	end
	if ((ischar(measureId) && (strcmpi("transfer", measureId) || strcmpi("all", measureId))) || \
	    ((measureId == 1) || (measureId == -1)))
		% Compute apparent transfer entropy
		if (measureParams.j == 0)
			error("Can't compute transfer entropy from a cell to itself (setting measureParams.j == 0)");
		end
		transferCalc = javaObject('infodynamics.measures.discrete.ApparentTransferEntropyCalculator', base, measureParams.k);
		transferCalc.initialise();
		transferCalc.addObservations(caStatesJInts, measureParams.j);
		avTransfer = transferCalc.computeAverageLocalOfObservations();
		printf("Average apparent transfer entropy (j=%d) = %.4f\n", measureParams.j, avTransfer);
		javaLocalValues = transferCalc.computeLocalFromPreviousObservations(caStatesJInts, measureParams.j);
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(javaLocalValues, options.plotOptions);
		if (options.saveImages)
			print(sprintf("figures/transfer-%d.eps", measureParams.j), "-color", "-deps");
		end
		plottedOne = true;
	end
	if ((ischar(measureId) && (strcmpi("transfercomplete", measureId) || strcmpi("completetransfer", measureId) || strcmpi("all", measureId))) || \
	    ((measureId == 2) || (measureId == -1)))
		% Compute complete transfer entropy
		if (measureParams.j == 0)
			error("Can't compute transfer entropy from a cell to itself (setting measureParams.j == 0)");
		end
		transferCalc = javaObject('infodynamics.measures.discrete.CompleteTransferEntropyCalculator', ...
			base, measureParams.k, neighbourhood - 2);
		transferCalc.initialise();
		% The offsets of the parents (see runCA for how this is computed, especially for even neighbourhood):
		fullSetOfParents = ceil(-neighbourhood / 2) : ceil(-neighbourhood / 2) + (neighbourhood-1);
		% Offsets of all parents can be included here - even 0 and j, these will be eliminated internally:
		transferCalc.addObservations(caStatesJInts, measureParams.j, octaveToJavaIntArray(fullSetOfParents));
		avTransfer = transferCalc.computeAverageLocalOfObservations();
		printf("Average complete transfer entropy (j=%d) = %.4f\n", measureParams.j, avTransfer);
		javaLocalValues = transferCalc.computeLocalFromPreviousObservations(caStatesJInts, ...
					measureParams.j, octaveToJavaIntArray(fullSetOfParents));
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(javaLocalValues, options.plotOptions);
		if (options.saveImages)
			print(sprintf("figures/transferComp-%d.eps", measureParams.j), "-color", "-deps");
		end
		plottedOne = true;
	end
	if ((ischar(measureId) && (strcmpi("separable", measureId) || strcmpi("all", measureId))) || \
	    ((measureId == 3) || (measureId == -1)))
		% Compute separable information
		separableCalc = javaObject('infodynamics.measures.discrete.SeparableInfoCalculator', ...
			base, measureParams.k, neighbourhood - 1);
		separableCalc.initialise();
		% The offsets of the parents (see runCA for how this is computed, especially for even neighbourhood):
		fullSetOfParents = ceil(-neighbourhood / 2) : ceil(-neighbourhood / 2) + (neighbourhood-1);
		% Offsets of all parents can be included here - even 0 and j, these will be eliminated internally:
		separableCalc.addObservations(caStatesJInts, octaveToJavaIntArray(fullSetOfParents));
		avSeparable = separableCalc.computeAverageLocalOfObservations();
		printf("Average separable information = %.4f\n", avSeparable);
		javaLocalValues = separableCalc.computeLocalFromPreviousObservations(caStatesJInts, ...
					octaveToJavaIntArray(fullSetOfParents));
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(javaLocalValues, options.plotOptions);
		if (options.saveImages)
			print(sprintf("figures/separable-k%d.eps", measureParams.k), "-color", "-deps");
		end
		plottedOne = true;
	end

	if (not(plottedOne))
		error(sprintf("Supplied measureId %s did not match any measurement types", measureId));
	end

	toc
end

