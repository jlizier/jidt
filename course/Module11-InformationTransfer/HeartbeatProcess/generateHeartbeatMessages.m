% function data = generateHeartbeatMessages(lambda1, lambda0, N)
%
% Generate sample time series data, of N time steps, for two processes.
% The first process s is a Poisson switching process, where the next value s_{n+1}
%  is a simple function of its previous value s_n:
%
%  p(s_{n+1} = 1 | s_n = 0) = lambda1
%  p(s_{n+1} = 0 | s_n = 0) = 1-lambda1
%  p(s_{n+1} = 1 | s_n = 1) = 1-lambda0
%  p(s_{n+1} = 0 | s_n = 1) = lambda0
%
% The second process t simply copies the previous value of s.
%
% We have for the steady-state probabilities:
% p(s_n = 0) = lambda0 / (lambda0 + lambda1)
% p(s_n = 1) = lambda1 / (lambda0 + lambda1)
% Proof:
% p(1) = p(0).p(s_{n+1} = 1 | s_n = 0) + p(1).p(s_{n+1} = 1 | s_n = 1)
%      = p(0).lambda1 + p(1).(1-lambda0)
% p(1) = lambda1/lambda0 . p(0)
% and: p(1) + p(0) = 1
% Solve these two equations to get the above result.
%
% Inputs:
% - lambda1 - probability of switching to 1 when current value is a 0
% - lambda0 - probability of switching to 0 when current value is a 1
% - N - number of time steps to generate
%
% Output:
% - data - Nx2 matrix. First column is process s, second column is process t.
%
% Copyright (C) 2017, Joseph T. Lizier
% Distributed under GNU General Public License v3
%
function data = generateHeartbeatMessages(lambda1, lambda0, N)

	s = zeros(N,1);
	t = zeros(N,1);
	% One can show:
	p1 = lambda1 ./ (lambda0 + lambda1);
	p0 = lambda0 ./ (lambda0 + lambda1);
	% Now assign s(1), t(1) according to these probabilities:
	s(1) = (rand() < p1) .* 1;
	t(1) = (rand() < p1) .* 1;
	% And assign the remaining time series activity based on the process
	%  dynamics.
	% The way I've done it here is rather slow, but fast enough for
	%  our purposes ... and v fast to write :)
	for n = 2:N
		% Assign next value conditioned on previous:
		if (s(n-1) == 0)
			% Prob that s(n) switches to 1 is lambda1
			s(n) = (rand() < lambda1) .* 1;
		else
			% Prob that s(n) stays at 1 is (1-lambda0)
			s(n) = (rand() < 1 - lambda0) .* 1;
		end
	end
	% Now do the delayed copy to the target
	t(2:N) = s(1:N-1);
	
	data(:,1) = s;
	data(:,2) = t;
end

