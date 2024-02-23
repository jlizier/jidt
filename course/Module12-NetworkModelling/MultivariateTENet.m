% Copyright (C) 2023, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week10/ca54.txt');

networkSize = 100;
network = zeros(networkSize);
% Compute for all targets:
for d = 1:networkSize
    
    fprintf('Beginning greedy selection of parents for %d\n', d);

    % Column indices start from 1 in Matlab:
    destination = octaveToJavaIntArray(data(:, d));

    conditionalSet = [];
    
    while (true)
    
        % Precondition: we have already selected the parents in
        % conditionalSet, now we check if we can add to this:
        
        results = zeros(1, networkSize);
        pValues = zeros(1, networkSize);
        
        % 1. Construct the calculator:
        if (length(conditionalSet) > 0)
            calc = javaObject('infodynamics.measures.discrete.ConditionalTransferEntropyCalculatorDiscrete', 2, 4, length(conditionalSet));
        else
            calc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 4, 1, 1, 1, 1);
        end
        % 2. No other properties to set for discrete calculators.

        for s = 1:networkSize
            % For each source-dest pair:
            if (s == d) || ~isempty(find(conditionalSet == s))
                pValues(s) = 1;
                continue;
            end
            % Column indices start from 1 in Matlab:
            source = octaveToJavaIntArray(data(:, s));
            conditional = octaveToJavaIntArray(data(:,conditionalSet));

            % 3. Initialise the calculator for (re-)use:
            calc.initialise();
            % 4. Supply the sample data:
            if (length(conditionalSet) > 0)
                calc.addObservations(source, destination, conditional);
            else
                calc.addObservations(source, destination);
            end
            % 5. Compute the estimate:
            result = calc.computeAverageLocalOfObservations();
            % 6. Compute the (statistical significance via) null distribution analytically:
            measDist = calc.computeSignificance();
            
            results(s) = result;
            pValues(s) = measDist.pValue;
        end
        
        % Check which was the strongest source:
        [maxTE, maxSourceIndex] = max(results);
        if (pValues(maxSourceIndex) < 0.05/(networkSize*(networkSize-1)))
            fprintf('Selected source %d, with TE(%d->%d | %s)=%.5f, p-value=%.5f (conditioning on %d parents)\n', ...
                maxSourceIndex, maxSourceIndex, d, strjoin(string(conditionalSet)), maxTE, ...
                pValues(maxSourceIndex), length(conditionalSet));
            % Add this new source:
            conditionalSet(end+1) = maxSourceIndex;
        else
            fprintf('-- Max TE source %d was not significant (TE(%d->%d | %s)=%.5f, p-value=%.6f), quitting\n', ...
                maxSourceIndex, maxSourceIndex, d, strjoin(string(conditionalSet)), maxTE, pValues(maxSourceIndex));
            break;
        end
    end
    
    % Postcondition: conditionalSet holds the parents for target d
    network(conditionalSet, d) = 1;
end

figure(); imagesc(network); xlabel('target'); ylabel('source'); title('Multivariate effective network via TE'); h = colorbar; ylabel(h, 'Connections');

