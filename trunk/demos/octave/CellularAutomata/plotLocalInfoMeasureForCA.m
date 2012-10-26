%
% 
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
%    - "active", 0 - active information storage (requires options.k)
%    - "all", -1 - plot all measures
% - measureParams - a structure containing options as described for each measure above:
%    - measureParams.k - history length for information dynamics measures
% - options - a stucture containing a range of other options, i.e.:
%    - plotOptions - structure as defined for the plotRawCa function
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
	javaaddpath('../../infodynamics.jar');

	% Simulate and plot the CA
	caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false);
	if (options.plotRawCa)
		plotRawCa(caStates, options.plotOptions, options.saveImages);
	end
	figNum = 2;
	toc

	% convert the states to a format usable by java:
	caStatesJInts = octaveToJavaIntMatrix(caStates);
	toc

	% Make the local information dynamics measurement
	if ((ischar(measureId) && (strcmp("active", measureId) || strcmp("all", measureId))) || \
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
	end

	toc
end

