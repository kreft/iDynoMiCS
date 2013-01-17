function plotAgentStructure(a,time,sName)
% plotAgentStructure(a,time,sName)
%
% this plots the agent structure given in 'a'
%
% to keep previous plot, include legend entry 'sName'
%

if nargin < 3
    sName = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first remove any species that have no individuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rem = [];
keep = [];
for i=1:numel(a)
	if size(a(i,1).data,1) > 0
		keep = [keep; i];
	else
		rem = [rem; i];
	end
end
for i=1:length(rem)
	fprintf('No individuals of species %s to plot.\n',a(rem(i)).name);
end
a = a(keep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nspec = numel(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort the species so that those with the most members are plotted first
%sorder = zeros(nspec,2);
%for i=1:nspec
%	sorder(i,1) = calcAgentData(a(i),'number');
%	sorder(i,2) = i;
%end
%sorder = sortrows(sorder,-1);
%anew = a;
%for i=1:nspec
%	anew(i) = a(sorder(i,2));
%end
%a = anew;


% create storage arrays for desired values
xyz = cell(nspec,1);
rad = cell(nspec,1);

% scale the radius values for more visible plots
radscale = 4;

% fill storage arrays
for i=1:nspec
	% copy position columns for all individuals
	xyz{i} = a(i).data(:,a(i).header.locationX:a(i).header.locationZ);

	% copy radius column for all individuals and scale the radius
	% (this if statement is used for backwards compatibility)
	if isfield(a(i).header,'totalRadius')
		rad{i} = a(i).data(:,a(i).header.totalRadius)*radscale;
	else
		rad{i} = a(i).data(:,a(i).header.radius)*radscale;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make agent names legible in plots
for i=1:nspec
    a(i).name = strrep(a(i).name,'_','\_');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(sName,'')
    hold on;
end

% this sets up the legend by plotting one cell of each type
% (put the points outside the domain)
for i=1:nspec
    plot(-10,-10,'o',...
		'MarkerSize',12,...
		'MarkerFaceColor',a(i).color,...
		'MarkerEdgeColor','k');
    hold on;
end
if strcmp(sName,'')
    legend(a.name);
else
    legend(sName,a.name);
end

% now go through and plot the points
maxht = 0;
maxrt = 0;
for i=1:nspec
    nparts = size(xyz{i},1);
    fprintf('Individuals of %s: %i\n',a(i).name,nparts);
    for j=1:nparts
        plot(xyz{i}(j,2),xyz{i}(j,1),'o',...
			'MarkerSize',rad{i}(j),...
			'MarkerFaceColor',a(i).color,...
			'MarkerEdgeColor','k');
        maxht = max(maxht,xyz{i}(j,1));
        maxrt = max(maxrt,xyz{i}(j,2));
    end
end

if strcmp(sName,'')
	axis equal;
	xlim([0 a(1).nJ*a(1).resolution])
	ylim([0 max([40,1.2*maxht])]);
else
	axis equal;
	xlim([0 a(1).nJ*a(1).resolution])
	ylim([0 a(1).nI*a(1).resolution])
end

xlabel('Y [\mum]');
ylabel('X [\mum]');
title(sprintf('Agents at %g Hours',time));

drawnow;

