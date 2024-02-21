% function stringRepresentation = translateMove(move)
%
% Returns a string representation of the given move index:
% 0 -> scissors
% 1 -> paper
% 2 -> rock

function stringRepresentation = translateMove(move)
	if (move == 0)
		stringRepresentation = 'scis';
	elseif (move == 1)
		stringRepresentation = 'papr';
	elseif (move == 2)
		stringRepresentation = 'rock';
	else
		error('Move %d not recognised', move);
	end
end

