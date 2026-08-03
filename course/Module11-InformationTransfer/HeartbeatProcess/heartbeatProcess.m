% Copyright (C) 2017-, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/infodynamics.jar');
% Add utilities to the path
addpath('/home/joseph/JIDT/infodynamics-dist-1.6.1/demos/octave');

% 0. Load/prepare the data:
data = generateHeartbeatMessages(0.05, 0.2, 100000);
% Column indices start from 1 in Matlab:
source = octaveToJavaIntArray(data(:,1));
destination = octaveToJavaIntArray(data(:,2));

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 1, 1, 1, 1, 1);
% 2. No other properties to set for discrete calculators.
% 3. Initialise the calculator for (re-)use:
calc.initialise();
% 4. Supply the sample data:
calc.addObservations(source, destination);
% 5. Compute the estimate:
result = calc.computeAverageLocalOfObservations();

fprintf('TE_Discrete(col_0 -> col_1) = %.4f bits\n', ...
	result);

locals = calc.computeLocalFromPreviousObservations(source, destination);
locals = javaMatrixToOctave(locals); % Conversion only required if running Octave (optional for Matlab)
figure(1); plot(data(1:100,1), 'rx'); hold on; plot(data(1:100,2), 'bo'); hold off; legend('source', 'target'); xlabel('n'); ylabel('state');
figure(2); plot(locals(1:100), 'rx'); ylabel('TE (bits)'); hold on; yyaxis right; plot(data(1:100,2), 'bo'); hold off; legend('locals', 'target'); xlabel('n'); ylabel('state');

%%
% Challenge task -- add lagged MI to the plot:

% Create MI calc with time lag 1
miCalc = javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2, 2, 1);
% 2. No other properties to set for discrete calculators.
% 3. Initialise the calculator for (re-)use:
miCalc.initialise();
% 4. Supply the sample data:
miCalc.addObservations(source, destination);
% 5. Compute the estimate:
result = miCalc.computeAverageLocalOfObservations();
localsMI = miCalc.computeLocalFromPreviousObservations(source, destination);
localsMI = javaMatrixToOctave(localsMI); % Conversion only required if running Octave (optional for Matlab)

% And our localsMI don't have a result for timestep 1 (where there is no history) yet TE does,
%  so align these properly:
localsMI = [0; localsMI];

% And plot them with the local TE values:
figure(3); plot(locals(1:100), 'rx'); ylabel('Information (bits)'); xlabel('n'); hold on;
plot(localsMI(1:100), 'g+');
yyaxis right; plot(data(1:100,2), 'bo'); hold off;
legend('local TE', 'local MI', 'target');
ylabel('state');
