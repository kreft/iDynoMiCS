function [adata,time] = getAgentData(nFile,type,skipeps)
% [adata,time] = getAgentData(nFile,type,skipeps)
%
% this is similar to the 'loadAgents' routine, but
% allows one to specify the data output kind
%
% nFile is the iterate number or 'last' for the last iterate
%
% type may be:
%	'State' for separate data for each agent
%
%	'Sum' for data summed by species type
%
%	'Grid' to return a solute-sized grid with the agent
%			information copied to it (total biomass per species)
%
%	'SplitGrid' to return a solute-sized grid with the agent
%			information split into active, inert, and EPS portions
%
% set skipeps to 1 to ignore any species with EPS or eps in their title

if nargin < 2
    type = 'State';
end

if nargin < 3
	skipeps = 0;
end

% default colors that are assigned to each species
acolor = {'r','b','g','y','m','c','o','k'};

% make lowercase for easier string comparisons
type = lower(type);

dogrid = 0;
dosplit = 0;
if strcmp(type,'grid') || strcmp(type,'splitgrid') || strcmp(type,'split')
    
    dogrid = 1;
    
    if strcmp(type,'splitgrid')
        dosplit = 1;
    end
    
    type = 'state';
end


if strcmp(type,'state') type='State'; end
if strcmp(type,'sum') type='Sum'; end


if isnumeric(nFile)
	stateFile = sprintf('agent_%s(%i).xml',type,nFile);
else
	stateFile = sprintf('agent_%s(%s).xml',type,nFile);
	if strcmp(nFile,'last')
		stateFile = sprintf('lastIter/agent_%s(last).xml',type);
	end
end

% decompress the file if needed
if ~exist(stateFile,'file')
    zipFile = sprintf('agent_%s.zip',type);
    success = unzipSelect(zipFile,stateFile);   

	if ~success
		fprintf('Could not read file "%s"\n',stateFile);
		fprintf('or unzip it from "%s".\n',zipFile);
		adata = 0;
		time = -1;
		return;
	end
end


% now starts the reading in of the XML data
xDoc = xmlread(stateFile);
xRoot = xDoc.getDocumentElement;

time = str2num(xRoot.getElementsByTagName('simulation').item(0).getAttribute('time'));

grid = xRoot.getElementsByTagName('grid').item(0);
if grid ~= []
	res = str2num(xRoot.getElementsByTagName('grid').item(0).getAttribute('resolution'));
	nI = str2num(xRoot.getElementsByTagName('grid').item(0).getAttribute('nI'));
	nJ = str2num(xRoot.getElementsByTagName('grid').item(0).getAttribute('nJ'));
	nK = str2num(xRoot.getElementsByTagName('grid').item(0).getAttribute('nK'));
else
	% since we don't have the resolution or gridpoint numbers, get them from the solute file
	[s] = getEnvData(nFile);
	res = s(1,1).resolution;
	nI = s(1,1).nI;
	nJ = s(1,1).nJ;
	nK = s(1,1).nK;
	clear s;
end

allSpecies = xDoc.getElementsByTagName('species');

adata = struct([]);

for iSpecies=1:allSpecies.getLength
    speciesName = char(allSpecies.item(iSpecies-1).getAttribute('name'));

	if skipeps && length(strfind(lower(speciesName),'eps')) > 0
		fprintf('Skipping species %s.\n',speciesName);
		continue;
	end

	% read the header but then create a record of which column holds
	% each data value
    header = char(allSpecies.item(iSpecies-1).getAttribute('header'));
	labels = regexp(header,'\w*','match');
	labelColumnMap = ([]);
	for i=1:length(labels)
		labelColumnMap = setfield(labelColumnMap,{1},labels{i},i);
	end

	% read in the data values
    x = allSpecies.item(iSpecies-1).getFirstChild;
    if ~isempty(x)
        x = str2num(x.getData);
    end

    adata = setfield(adata,{iSpecies,1},'name',speciesName);
    adata = setfield(adata,{iSpecies,1},'color',acolor{iSpecies});
    adata = setfield(adata,{iSpecies,1},'resolution',res);
    adata = setfield(adata,{iSpecies,1},'nI',nI);
    adata = setfield(adata,{iSpecies,1},'nJ',nJ);
    adata = setfield(adata,{iSpecies,1},'nK',nK);
    adata = setfield(adata,{iSpecies,1},'header',labelColumnMap);
    adata = setfield(adata,{iSpecies,1},'data',x);
end

% when we want the data on a solute-style grid
if dogrid
    adata = agentListToGrid(adata,dosplit);
end

