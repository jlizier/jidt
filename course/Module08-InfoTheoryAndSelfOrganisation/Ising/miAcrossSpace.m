% This script examines how MI between cells changes with spatial distance
%  between the cells (explored via distance across the x-axis) -- at which
%  temperatures does the system exhibit long-correlations?
% The script also examines how this changes with temperature, to explore
%  the phase transition as well. Can you see the effect of the phase
%  transition brought out here? (For finite-sized systems, the transition
%  point is usually offset from the theoretical critical point to higher
%  temperatures).

% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% We may want to loop over the data files later:
tempsToExamine = [2.00, 2.20, 2.26, 2.27, 2.28, 2.34, 2.54];
labels = {};

% Check one data file - the way we reshape it into a grid may be useful for you to copy later:
data = load('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ising/spins-2.27.txt');
modelSize = sqrt(size(data,2));
timeSteps = size(data, 1);
dataInLayout = zeros(timeSteps, modelSize, modelSize);
for t = 1:timeSteps
    dataInLayout(t,:,:) = reshape(data(t,:), modelSize, modelSize);
end
% Let's plot the last time step to check that it looks ok:
imagesc(squeeze(dataInLayout(timeSteps,:,:)));
title(sprintf('Spins at t=%d for T=2.27', timeSteps));

% Now compute and plot MI vs x-separation for each temp:
figure();
hold off;
for ti = 1:length(tempsToExamine)
    labels{ti} = sprintf('T=%.2f', tempsToExamine(ti));

    % 0. Load/prepare the data:
    data = load(sprintf('/home/joseph/TeachingPlayground/CSYS5030/Week8/Ising/spins-%.2f.txt', tempsToExamine(ti)));
    
    % Let's reshape the data so it's easier to work with
    dataInLayout = zeros(timeSteps, modelSize, modelSize);
    for t = 1:timeSteps
        dataInLayout(t,:,:) = reshape(data(t,:), modelSize, modelSize);
    end
    
    misVsXDist = zeros(1, ceil(modelSize/2));
    % 1. Construct the calculator:
    calc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2, 2, 0);
    % 2. No other properties to set for discrete calculators.
    for xDiff = 1:ceil(modelSize/2)
        % We'll calculate MI for cells separated by xDiff
        % 3. Initialise the calculator for (re-)use:
        calc.initialise();
        for x = 1:modelSize
            for y = 1:modelSize
                % Column indices start from 1 in Matlab:
                source = dataInLayout(:,x,y);
                destination = dataInLayout(:,mod(x+xDiff - 1, modelSize)+1,y);
                % 4. Supply the sample data:
                calc.addObservations(source, destination);
            end
        end
        % 5. Compute the estimate:
        misVsXDist(xDiff) = calc.computeAverageLocalOfObservations();
        % fprintf('MI_Discrete(temp=%.2f,xDiff=%d) = %.4f bits\n', ...
	    %    tempsToExamine(ti), xDiff, misVsXDist(xDiff));
    end
    
    plot(misVsXDist, 'x')
    hold on;
end
hold off;
legend(labels);
title('MI between cells as a function of x-separation');
ylabel('MI (bits)');
xlabel('x difference (cells)')
