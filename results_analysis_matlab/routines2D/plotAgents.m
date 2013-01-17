function plotAgents(nFile,sName)
% plotAgents(nFile,sName)
%
% input sName (solute name) to retain previous plot
% (assumes the solute was plotted previously)
%

if nargin < 1
	help plotAgents;
	return;
end

if nargin < 2
    sName = '';
else
	% make solute name legible in plots
	sName = strrep(sName,'_','\_');
end

[a,time] = getAgentData(nFile,'State');

% now just call the generic agent plotting routine
plotAgentStructure(a,time,sName);

