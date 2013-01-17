function plotContour(nFile,fName)
% plotContour(nFile,fName)
%
% make a contour plot of data field fName at iterate nFile
%
% run showRunInfo to see possible fields to plot.
% you may also use 'allbio' to plot a field for total biomass,
% and 'interface' to plot the biofilm/liquid interface
%

if nargin < 1
	help plotContour;
	return;
end

fs = 24;
lw = 3;


[s,a,r,time] = getEnvData(nFile);
nsols = numel(s);
nspec = numel(a);
nreac = numel(r);

a = getAgentData(nFile,'Grid');
nspec = size(a,1);
ntype = size(a,2);

if ~exist('fName','var')
	fName = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we need to find out what we're plotting, so check each in turn
%
% these routines just march through each field and check whether
% the input name matches any of the field names, and if it does
% then it copies the needed data out for plotting

thetype = '';

is = 1;
while is <= nsols && ~strcmp(s(is).name,fName)
    is = is+1;
end
if is <= nsols
    thetype = 'solute';
    fUnits = s(is).unit;
    res = s(is).resolution;
    D = s(is).data;
end

if strcmp(thetype,'')
    ntype = size(a,2);
    is = 1;
    it = 1;
    while is <= nspec && ~strcmp(a(is,it).name,fName)
        it = it+1;
        if it > ntype
            is = is+1;
            it = 1;
        end
    end
    if is <= nspec
        thetype = 'agent';
        fUnits = a(is,it).unit;
        res = a(is,it).resolution;
        D = a(is,it).data;
    end
end

if strcmp(thetype,'')
    is = 1;
    while is <= nreac && ~strcmp(r(is).name,fName)
        is = is+1;
    end
    if is <= nreac
        thetype = 'reaction';
        fUnits = r(is).unit;
        res = r(is).resolution;
        D = r(is).data;
    end
end

% add ability to plot ALL biomass
if strcmp(thetype,'')
    if strcmp(fName,'allbio')
        fName = 'Total Biomass';
        thetype = 'agent';
        fUnits = a(1).unit;
        res = a(1).resolution;
        D = calcBiofilmData(a,'biomass');
    end
end

% plot an interface contour line
if strcmp(thetype,'')
    if strcmp(fName,'interface')
        fName = 'Interface';
        thetype = 'interface';
        fUnits = '-';
        res = a(1).resolution;
        D = calcBiofilmData(a,'biomass');
		D = ones(size(D)).*(D > 0);
    end
end

if strcmp(thetype,'')
    % we couldn't find the name
    fprintf('Field %s does not exist.\n',fName);
    fprintf('Choose from:\n');
    showRunInfo(nFile);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make reaction name legible in plots
fName = strrep(fName,'_','\_');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the contour plot
figure;

x = ((1:size(D,2)) - .5)*res;
y = ((1:size(D,1)) - .5)*res;
[X,Y] = meshgrid(x,y);

if strcmp(thetype,'interface')
	contour(X,Y,D,[0 0],'-k','LineWidth',lw);
else
	contourf(X,Y,D);
	colorbar;
end
axis equal;
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
xlabel('Y [\mum]','FontSize',fs);
ylabel('X [\mum]','FontSize',fs);
title(sprintf('%s [%s]',fName,fUnits),'FontSize',fs);
set(gca,'FontSize',fs);


