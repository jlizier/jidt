% J. Lizier, 2012.
%
% Function to create a colourmap of a given length.
% 
% The largest values are given the darkest shades of only the primary colour (blue or red),
%  with fracPrim scaling the amount of the vector this occupies.
% The next values are given lighter shades of the primary colour, using green to
%  make it lighter; fracWithSecondary scales how much of the vector this occupies.
% The last values are given the lightest (closest to white) shades, using the last
%  colour to do this.
% 
% Inputs:
% - vLength - length of the colourmap
% - primaryIsBlue - whether we're making a blue or red colourmap
% - fracPrim - what proportion to make the darkest shade of the primary colour (default .25)
% - fracSecondary - what proportion to make the lighter shades with green added (default .4)
% - gammaPower - scale the colours non-linearly with this exponent. I think this might be a gamma correction but am unsure.
%
% vLength is the length of the colormap
function rgbmap = prepareColourmap(vLength, primaryIsBlue, fracPrim, fracWithSecondary, gammaPower)

	if (nargin < 5)
		gammaPower = 1;
	end
	if (nargin < 4)
		fracWithSecondary = .4;
	end
	if (nargin < 3)
		fracPrim = .25;
	end
	if (nargin < 2)
		primaryIsBlue = true;
	end
	if (nargin < 1)
		vLength = 64;
	end

	% Set the colormap to have blue for positive, scaled to the max of our local values
	% If the user wants red, we'll switch it around at the last minute.
	blueOnlyLength = floor(vLength .* fracPrim);
	greenAndBlueLength = floor(vLength .* fracWithSecondary);
	redOnLength = vLength - blueOnlyLength - greenAndBlueLength;
	bluevector = ((2.*blueOnlyLength:-1:blueOnlyLength+1)' ./ 2 ./ blueOnlyLength).^gammaPower; % Countdown from (1 to 0.5).^gammaPower
	bluevector = [ones(vLength - blueOnlyLength, 1); bluevector];
	greenvector = ((greenAndBlueLength:-1:1)' ./ greenAndBlueLength).^gammaPower; % Count down from (1 to 0).^gammaPower
	greenvector = [ones(redOnLength, 1); greenvector; zeros(blueOnlyLength, 1)];
	redvector = ((redOnLength:-1:1)' ./ redOnLength).^gammaPower; % Count down from (1 to 0).^gammaPower
	redvector = [redvector; zeros(vLength - redOnLength, 1)];

	if (primaryIsBlue)
		rgbmap = [ redvector, greenvector, bluevector];
	else
		% Primary colour is red, so swap the blue and red columns here
		rgbmap = [ bluevector, greenvector, redvector];
	end
end

