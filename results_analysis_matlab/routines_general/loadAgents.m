function data=loadAgents(i,deleteAfterUnZip)

error(nargchk(1,2,nargin));
if(nargin==1)
    deleteAfterUnZip=false;
end

if(i==-1) index='last'; else index = num2str(i);end
fileName = ['agent_State(',num2str(index),').xml'];

test=exist(fileName,'file')~=0;
flagUnzip = false;

% Try to decompress the file if it possible
if ~test & exist('agent_State.zip','file')~=0
    flagUnzip=unzipSelect('agent_State.zip',fileName);   
    if flagUnzip
        disp(['File ',fileName,' decompressed']);
    else
        error('File not found');
    end
end

xDoc = xmlread(fileName);

if deleteAfterUnZip&flagUnzip
    delete(fileName);
    disp(['File ',fileName,' deleted']);
end

xRoot = xDoc.getDocumentElement;


time = xRoot.getElementsByTagName('simulation').item(0).getAttribute('time');

allSpecies = xDoc.getElementsByTagName('species');
data = struct([]);

for iSpecies=1:allSpecies.getLength
    speciesName = char(allSpecies.item(iSpecies-1).getAttribute('name'));
    header = char(allSpecies.item(iSpecies-1).getAttribute('header'));
    x = allSpecies.item(iSpecies-1).getFirstChild;
    if ~isempty(x)
        x = str2num(x.getData);
    end
    
    data = setfield(data,{iSpecies},'name',speciesName);
    data = setfield(data,{iSpecies},'header',header);
    data = setfield(data,{iSpecies},'data',x);    
end