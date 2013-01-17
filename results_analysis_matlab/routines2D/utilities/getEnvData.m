function [s,a,r,time,heights,p] = getEnvData(nFile)
% [sols,agents,reacs,time,heights,pressure] = getEnvData(nFile)
%
% load environment data for the given iterate
%
% heights is an array with the mean, stddev, and max interface height data
%
%


if isnumeric(nFile)
	stateFile = sprintf('env_State(%i).xml',nFile);
else
	stateFile = sprintf('env_State(%s).xml',nFile);
	if strcmp(nFile,'last')
		stateFile = 'lastIter/env_State(last).xml';
	end
end

% decompress the file if needed
if ~exist(stateFile,'file')
    zipFile = 'env_State.zip';
    success = unzipSelect(zipFile,stateFile);

	if ~success
		fprintf('Could not read or unzip file "%s".\n',stateFile);
		s = 0;
		a = 0;
		r = 0;
		time = -1;
		heights = 0;
		return;
	end
end


% now starts the reading in of the XML data
xDoc = xmlread(stateFile);
xRoot = xDoc.getDocumentElement;
data = struct([]);
time = str2num(xRoot.getElementsByTagName('simulation').item(0).getAttribute('time'));

% read in the height information
heights = -ones(3,1);
thickRoot = xRoot.getElementsByTagName('thickness').item(0);
if thickRoot ~= []
	heights(1) = str2num(thickRoot.getElementsByTagName('mean').item(0).getFirstChild.getData);
	heights(2) = str2num(thickRoot.getElementsByTagName('stddev').item(0).getFirstChild.getData);
	heights(3) = str2num(thickRoot.getElementsByTagName('max').item(0).getFirstChild.getData);
end


% this is all copied from the 'loadSolute' routine from Laurent,
% as is the 'isShaped' subroutine below

allSolutes = xDoc.getElementsByTagName('solute');

for iSolute=1:allSolutes.getLength
    soluteName = char(allSolutes.item(iSolute-1).getAttribute('name'));
    soluteUnit = char(allSolutes.item(iSolute-1).getAttribute('unit'));
	if length(soluteUnit)==0
		% this is for the old-style unit output, kept for compatibility
		soluteUnit = soluteName(findstr('(',soluteName)+1:findstr(')',soluteName)-1);
		soluteName = soluteName(1:findstr('(',soluteName)-1);
	end

    nI = str2num(allSolutes.item(iSolute-1).getAttribute('nI'));
    nJ = str2num(allSolutes.item(iSolute-1).getAttribute('nJ'));
    nK = str2num(allSolutes.item(iSolute-1).getAttribute('nK'));
    res = str2num(allSolutes.item(iSolute-1).getAttribute('resolution'));
    x=str2num(allSolutes.item(iSolute-1).getFirstChild.getData);
    

    data = setfield(data,{iSolute},'name',soluteName);
    data = setfield(data,{iSolute},'unit',soluteUnit);
    data = setfield(data,{iSolute},'resolution',res);
    data = setfield(data,{iSolute},'nI',nI);
    data = setfield(data,{iSolute},'nJ',nJ);
    data = setfield(data,{iSolute},'nK',nK);
    data = setfield(data,{iSolute},'data',isShaped(x,nI,nJ,nK));
    %data = setfield(data,{iSolute},'data', x);
end

% now to split the data structure into structures with solutes (s),
% agent data (a) and reaction data (r), we must do some trickery
%
% go through the solutes structure to find the agent info, assuming the following:
% 1. the solutes are the first nsols entries in the structure
% 2. all reactions have '-rate' appended to their name
% 3. the biomass entries appear after the ractions and do NOT have '-rate'
% 4. the last entry is for the totalBiomass value

% this finds number of solutes
is=1;
while length(strfind(data(is).name,'-rate')) == 0
    is = is+1;
end
nsols = is-1;
s = data(1:nsols);

% this finds the start of the agent data by looking for where
% '-rate' is not part of the name
is=nsols+1;
while length(strfind(data(is).name,'-rate')) > 0
    is = is+1;
end

% now we can copy the agent data out of the solute structure
a = data(is:numel(data)-1);
%nspec = numel(a);

% and we can copy the reaction data out too
r = data(nsols+1:is-1);
% also add time element to the reaction units if it's not there
% THIS ASSUMES A 1/hr UNIT OF TIME!!!
for ir=1:numel(r)
	if length(strfind(r(ir).unit,'hr-1'))==0
		r(ir).unit = sprintf('%s.hr-1',r(ir).unit);
	end
end
nreac = numel(r);

% if the pressure is included in the solutes, copy it to its own list
for is=1:nsols
	if length(strfind(data(is).name,'pressure')) > 0
		p = s(is);
		s(is) = [];
		break;
	end
end

% finally transpose the data to make it in columnar form
s = s';
a = a';
r = r';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%_______________________ Subfunction ____________________________
function Xnew=isShaped(x,nI,nJ,nK)

if(length(x))==(nI*nJ)
    % No padding
	% The dimension order is nJ,nI,nK because the data are written with
	% the k-dimension along the columns, and with the j-dimension varying
	% faster down the rows than the i-dimension. Values are placed by
	% column by the reshape() routine, so it will fill the j-dimension first
	% and i-dimension second. We then permute to shape it in regular I,J,K
	% order.
	x = reshape(x,nJ,nI);
 	Xnew = permute(x,[2,1]);
elseif (length(x))==(nI+2)*(nJ+2)
    %padding
	% The dimension order is nJ,nI,nK because the data are written with
	% the k-dimension along the columns, and with the j-dimension varying
	% faster down the rows than the i-dimension. Values are placed by
	% column by the reshape() routine, so it will fill the j-dimension first
	% and i-dimension second. We then permute to shape it in regular I,J,K
	% order.
	x = reshape(x,nJ+2,nI+2);
	x = permute(x,[2,1]);
    Xnew = x(2:end-1,2:end-1);
else    
    Xnew = x;
end
