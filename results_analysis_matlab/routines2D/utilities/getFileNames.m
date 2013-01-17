function [files] = getFileNames(thedir,inctype,exctype)
% [files] = getFileNames(thedir,inctype,exctype)
%
% return cell array of file names, with the option to only
% return those that have 'inctype' in their name and exclude
% any with 'exctype' in their name
%

ofiles = dir(thedir);
files = {};

% ignore '.' and '..' that are always the first two entries
for i=3:length(ofiles)
	% only for files
	if ~ofiles(i).isdir
		% only include directories that match inctype and
		% don't have exctype
		if nargin >= 2
			if strfind(ofiles(i).name,inctype)
				if nargin == 3
					if length(strfind(ofiles(i).name,exctype)) == 0
						files{length(files)+1,1} = ofiles(i).name;
					end
				else
					files{length(files)+1,1} = ofiles(i).name;
				end
			end
		else
			files{length(files)+1,1} = ofiles(i).name;
		end
	end
end

