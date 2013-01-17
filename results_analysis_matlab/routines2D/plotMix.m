function plotMix(nFile,sName)
% plotMix(nFile,sName)
%
% plots agents and the solute sName for the given iterate nFile
% (the options for sName are the same as in plotContour)
%

if nargin < 1
	help plotMix;
	return;
end

if nargin < 2
	fprintf('Need to specify contour to plot.\n');
	showRunInfo(nFile);
	return;
end

plotContour(nFile,sName);
plotAgents(nFile,sName);
