function showRunInfo(nFile)
% showRunInfo(nFile)
%
% print out data for the run at iterate nFile
%
% (if nFile is not given, the last iterate is used instead)
%

import java.util.zip.ZipFile;

if nargin < 1
    nFile = 'last';
end

[s,a,r,time] = getEnvData(nFile);
a = getAgentData(nFile,'Sum');
b = getBulkData(nFile);

fprintf('Time: %g\n',time);

fprintf('Bulks:\n');
for is=1:b(1,1).number_bulks
    fprintf('\t%s\n',b(is,1).bulk_name);
end

fprintf('Solutes:\n');
for is=1:length(s)
    fprintf('\t%s\n',s(is).name);
end

fprintf('Agents:\n');
for is=1:size(a,1)
    for js=1:size(a,2)
        fprintf('\t%s\n',a(is,js).name);
    end
end

fprintf('Reactions:\n');
for is=1:length(r)
    fprintf('\t%s\n',r(is).name);
end

% show list of iterates saved in the zip files
itlist = [];
zipfile = ZipFile('env_Sum.zip');
enumeration = zipfile.entries;
while enumeration.hasMoreElements
	xmlfile = enumeration.nextElement.toString.toCharArray';
	itlist = [itlist, str2num(xmlfile([9:strfind(xmlfile,')')-1]))];
end
zipfile.close;

itlist = sort(itlist);
fprintf('Available iterates:\n');
for i=1:length(itlist)
	fprintf('%i ',itlist(i));
	if mod(i,15)==0
		fprintf('\n');
	end
end
fprintf('\n');

