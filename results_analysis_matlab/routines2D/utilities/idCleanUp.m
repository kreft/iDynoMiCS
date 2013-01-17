function idCleanUp(dirs)
% takes a cell array of directories and cleans out the following files:
%
% pass in 'sub' to get all the subdirectories of this one
% pass in 'recursive' to clean ALL subdirectories recursively
%
%  env_Sum(*).xml
%  env_State(*).xml
%  agent_Sum(*).xml
%  agent_State(*).xml
%  it(*).pov
%

if nargin < 1
	dothecleaning;
	return;
end

if strcmp(dirs,'sub')
	dirs = getDirNames('.');
end

recursive = 0;
if strcmp(dirs,'recursive')
	recursive = 1;
	dirs = getDirNames('.');
end

for id=1:numel(dirs)
	if strcmp(dirs{id},'lastIter')
		%fprintf('Skipping lastIter\n');
		continue;
	end

	cd(dirs{id});

	fprintf('Cleaning %s [%i/%i]\n',dirs{id},id,numel(dirs));
	dothecleaning;

	if recursive
		idCleanUp('recursive');
	end

	cd ..;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dothecleaning();

delete('env_Sum(*).xml');
delete('env_State(*).xml');
delete('agent_Sum(*).xml');
delete('agent_State(*).xml');
delete('it(*).pov');


