function iters = getListOfIterates(fkind,fext)
% iters = getListOfIterates(fkind,fext)
%
% this routine gets a list of the iterates in the current directory
%
% if fext is not given, it assumes 'xml'
%
% (assumes names of the form 'agent_State(123).xml' and so on)

if nargin < 2
    fext = 'xml';
end

files = dir(sprintf('%s(*).%s',fkind,fext));
list = {};

% now we want to extract the numbers from the list of output files
% assume a name of the form: 'fkind(123).fext)'; we want 123 from that string
istart = length(fkind)+2;
for i=1:length(files)
	iend = strfind(files(i).name,')')-1;
	thenum = files(i).name(istart:iend);
	if ~strcmp(thenum,'last')
		% only add the number if it isn't a 'last' iterate
		list{end+1,1} = thenum;
	end
end

% finally prepare them for output
iters = sort(str2num(char(list)));

