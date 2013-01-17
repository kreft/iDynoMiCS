function [data,time,heights] = getBulkData(nFile)
% [b,time,heights] = getBulkData(nFile)
%
% load bulk data from env_Sum for the given iterate
%
% heights is an array with the mean, stddev, and max interface height data
%
% (note that the pressure field is removed by this routine; look at the code
%  for the relevant section that controls this at the end of the routine)
%

if isnumeric(nFile)
	stateFile = sprintf('env_Sum(%i).xml',nFile);
else
	stateFile = sprintf('env_Sum(%s).xml',nFile);
	if strcmp(nFile,'last')
		stateFile = 'lastIter/env_Sum(last).xml';
	end
end

% decompress the file if needed
if ~exist(stateFile,'file')
    zipFile = 'env_Sum.zip';
    success = unzipSelect(zipFile,stateFile);

	if ~success
		fprintf('Could not read or unzip file "%s".\n',stateFile);
		data = 0;
		time = -1;
		heights = 0;
		return;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%

% now starts the reading in of the XML data
xDoc = xmlread(stateFile);
xRoot = xDoc.getDocumentElement;
time = str2num(xRoot.getElementsByTagName('simulation').item(0).getAttribute('time'));

% read in the height information
heights = -ones(3,1);
thickRoot = xRoot.getElementsByTagName('thickness').item(0);
if thickRoot ~= []
	heights(1) = str2num(thickRoot.getElementsByTagName('mean').item(0).getFirstChild.getData);
	heights(2) = str2num(thickRoot.getElementsByTagName('stddev').item(0).getFirstChild.getData);
	heights(3) = str2num(thickRoot.getElementsByTagName('max').item(0).getFirstChild.getData);
end


%%%%%%%%%%%%%%%%%%%%%%%%% BULK DATA %%%%%%%%%%%%%%%%%%%%%%%%%

data = struct([]);
% this was created by BRIM to allow looking at bulk values

allBulks = xDoc.getElementsByTagName('bulk');

% this flag will be set if there is a pressure field
haspressure = 0;

for iBulks=1:allBulks.getLength
    bulkName = char(allBulks.item(iBulks-1).getAttribute('name'));

    allSolutes = allBulks.item(iBulks-1).getElementsByTagName('solute');
    allRates = allBulks.item(iBulks-1).getElementsByTagName('uptake_rate');

    for iSolutes=1:allSolutes.getLength
        data = setfield(data,{iBulks,iSolutes},'number_bulks',allBulks.getLength);
        data = setfield(data,{iBulks,iSolutes},'number_solutes',allSolutes.getLength);

        soluteName = char(allSolutes.item(iSolutes-1).getAttribute('name'));
        soluteUnit = char(allSolutes.item(iSolutes-1).getAttribute('unit'));

		if strcmp(soluteName,'pressure')
			haspressure = 1;
		end

        data = setfield(data,{iBulks,iSolutes},'bulk_name',bulkName);
        data = setfield(data,{iBulks,iSolutes},'solute_name',soluteName);
        data = setfield(data,{iBulks,iSolutes},'solute_unit',soluteUnit);

        x = allSolutes.item(iSolutes-1).getFirstChild;
        if ~isempty(x)
            x = str2num(x.getData);
        end
        data = setfield(data,{iBulks,iSolutes},'solute_data',x);

        
        % to get the rate data, use the same solute index but copy data
        % from the rates structure instead
        uptakeUnit = char(allRates.item(iSolutes-1).getAttribute('unit'));
        data = setfield(data,{iBulks,iSolutes},'uptake_unit',uptakeUnit);
        
        x = allRates.item(iSolutes-1).getFirstChild;
        if ~isempty(x)
            x = str2num(x.getData);
        end
        data = setfield(data,{iBulks,iSolutes},'uptake_data',x);
   end
end

% now remove the pressure
% (comment this section out if you want to keep the pressure field in)
if haspressure
	for iBulks=1:allBulks.getLength
		for iSolutes=1:allSolutes.getLength
			data = setfield(data,{iBulks,iSolutes},'number_solutes',allSolutes.getLength - 1);
			if strcmp(data(iBulks,iSolutes).solute_name,'pressure')
				data(:,iSolutes) = [];
			end
		end
	end
end


