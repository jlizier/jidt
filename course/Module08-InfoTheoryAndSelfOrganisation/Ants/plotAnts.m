
positionsx = load('/home/joseph/Work/Teaching/USyd/CSYS5030-InfoTheoryAndSelfOrg/Lectures/Module6-InfoTheoryAndSelfOrganisation/positionsx.txt');
positionsy = load('/home/joseph/Work/Teaching/USyd/CSYS5030-InfoTheoryAndSelfOrg/Lectures/Module6-InfoTheoryAndSelfOrganisation/positionsy.txt');
headings = load('/home/joseph/Work/Teaching/USyd/CSYS5030-InfoTheoryAndSelfOrg/Lectures/Module6-InfoTheoryAndSelfOrganisation/headings.txt');

figure();
for t = 1 :size(positionsx, 1)
	scatter(positionsx(t,:), positionsy(t,:), 'x');
	axis([-40, 40, -40, 40])
	title(sprintf('Time %d', t));
	pause(0.1)
end
