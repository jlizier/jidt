% function stringRepresentation = translateResult(gameResult)
%
% Returns a string representation of the given game result:
% -1 -> lose
% 0 -> tie
% 1 -> win

function stringRepresentation = translateResult(gameResult)
	if (gameResult == -1)
		stringRepresentation = 'los';
	elseif (gameResult == 0)
		stringRepresentation = 'tie';
	elseif (gameResult == 1)
		stringRepresentation = 'win';
	else
		error('Result %d not recognised', gameResult);
	end
end
