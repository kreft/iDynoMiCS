function [dirs,updir] = getDirNames(thedir,inctext,exctext)
% [dirs,updir] = getDirNames(thedir,inctext,exctext)
%
% return cell array of directory names, with the option to only
% return those that have 'inctext' in their name and do not include
% 'exctext' in their name
%
% also returns 'updir', which is a string like '../..' to return
% from one of the 'dirs' returned (uses first in the list as example)
%

if nargin < 1
	thedir = '.';
end

odirs = dir(thedir);
dirs = {};

for i=1:length(odirs)
	% ignore '.' and '..'
	if strcmp(odirs(i).name,'.') || strcmp(odirs(i).name,'..')
		continue;
	end

	% only for directories
	if odirs(i).isdir
		% only include directories that match inctext and
		% don't have exctext
		if nargin >= 2
			if strfind(odirs(i).name,inctext)
				if nargin == 3
					if length(strfind(odirs(i).name,exctext)) == 0
						dirs{length(dirs)+1,1} = odirs(i).name;
					end
				else
					dirs{length(dirs)+1,1} = odirs(i).name;
				end
			end
		else
			dirs{length(dirs)+1,1} = odirs(i).name;
		end
	end
end

% now create updir by looking at how many subdirectories are in the path
updir = '..';
if length(dirs) > 1
	nup = length(strfind(dirs{1},'/'));
	for i=1:nup
		updir = sprintf('%s/..',updir);
	end
end

