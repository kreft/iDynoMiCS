function data = loadSolute(i,deleteAfterUnZip)
% Project iDynoMiCS (For copyright information, please consult copyright.txt)
% Open a resultfile describing fields of concentration
%______________________________________________________
% Author : Laurent Lardon

error(nargchk(1,2,nargin));
if(nargin==1)
    deleteAfterUnZip=false;
end

if(i==-1) index='last'; else index = num2str(i);end
fileName = ['env_State(',num2str(index),').xml'];

test=exist(fileName,'file')~=0;
flagUnzip = false;

% Try to decompress the file if it is possible
if ~test & exist('env_State.zip','file')~=0
    flagUnzip=unzipSelect('env_State.zip',fileName);   
    if flagUnzip
        disp(['File ',fileName,' decompressed']);
    else
        error('File not found');
    end
end


xDoc = xmlread(fileName);
xRoot = xDoc.getDocumentElement;
time = xRoot.getElementsByTagName('simulation').item(0).getAttribute('time');

allSolute = xDoc.getElementsByTagName('solute');
data = struct([]);
for iSolute=1:allSolute.getLength
    soluteName = char(allSolute.item(iSolute-1).getAttribute('name'));
    nI = str2num(allSolute.item(iSolute-1).getAttribute('nI'));
    nJ = str2num(allSolute.item(iSolute-1).getAttribute('nJ'));
    nK = str2num(allSolute.item(iSolute-1).getAttribute('nK'));
    res = str2num(allSolute.item(iSolute-1).getAttribute('resolution'));
    x=str2num(allSolute.item(iSolute-1).getFirstChild.getData);

    data = setfield(data,{iSolute},'name',soluteName);
    data = setfield(data,{iSolute},'resolution',res);
    data = setfield(data,{iSolute},'data',isShaped(x,nI,nJ,nK));
end


%_______________________ Subfunction ____________________________
function X=isShaped(x,nI,nJ,nK)
nX=size(x,1);
nY=size(x,2);
nZ=size(x,3);

if(nX*nY*nZ)==(nI*nJ*nK)
    % No padding
	x = reshape(x,nJ,nI,nK);
	X = permute(x,[2,1,3]);
elseif (nX*nY*nZ)==(nI+2)*(nJ+2)*(nK+2)
    %padding
	x = reshape(x,nJ+2,nI+2,nK+2);
	x = permute(x,[2,1,3]);
    X = x(2:end-1,2:end-1,2:end-1);
else    
    X   = x;
end
