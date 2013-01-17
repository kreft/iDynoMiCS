function agrid = agentListToGrid(adata,dosplit)
% agrid = agentListToGrid(adata,dosplit)
%
% this takes individual agent data and converts it to a grid
%
% set dosplit to 1 to split active, inert, and EPS portions
%

if nargin < 2
    dosplit = 0;
end

nspec = numel(adata);
ntype = 1;

% get voxel volume (cube of the resolution)
res = adata(1,1).resolution;
vol = res^3;

% now copy new agent data in & create new 'agrid' fields
agrid = ([]);
for ispec=1:nspec
    agrid(ispec,1).name       = adata(ispec).name;
    agrid(ispec,1).resolution = adata(1,1).resolution;
    agrid(ispec,1).nI         = adata(1,1).nI;
    agrid(ispec,1).nJ         = adata(1,1).nJ;
    agrid(ispec,1).nK         = adata(1,1).nK;
    agrid(ispec,1).unit       = 'g.L-1'; % THIS IS ASSUMED!!
    agrid(ispec,1).data       = zeros(adata(1,1).nI,adata(1,1).nJ,adata(1,1).nK);
end

if dosplit
    for ispec=1:nspec
		% iterate over each species and copy the data in each field so that
		% we have it in triplicate, with the spots to be replaced by active
		% biomass, inert biomass, and EPS
		agrid(ispec,2) = agrid(ispec,1);
		agrid(ispec,3) = agrid(ispec,1);

		% now modify the name of each component
        specname = adata(ispec).name;
        agrid(ispec,1).name = sprintf('%s-biomass',specname);
        agrid(ispec,2).name = sprintf('%s-inert',specname);
        agrid(ispec,3).name = sprintf('%s-capsule',specname);
    end
end

% now that we have the spots ready, iterate over agents and place the
% data into the correct field data
for ispec=1:nspec

	% data locations for this species
	if isfield(adata(ispec).header,'locationX')
		locx = adata(ispec).header.locationX;
		locy = adata(ispec).header.locationY;
		locz = adata(ispec).header.locationZ;
	else
		fprintf('*** Using zero location for species %s. ***\n',adata(ispec).name);
		locx = 0;
		locy = 0;
		locz = 0;
	end

	havebio = isfield(adata(ispec).header,'biomass');
	if havebio
		locb = adata(ispec).header.biomass;
	end

	havenrt = isfield(adata(ispec).header,'inert');
	if havenrt
		loci = adata(ispec).header.inert;
	end

	havecap = isfield(adata(ispec).header,'capsule');
	if havecap
		locc = adata(ispec).header.capsule;
	end

	% now iterate over all the agents in the structure and copy them to the grid
    for iagent=1:size(adata(ispec).data,1)
        drow = adata(ispec).data(iagent,:);

        % first get the location on the grid
        i = floor(drow(locx)/res)+1;
        j = floor(drow(locy)/res)+1;
        k = floor(drow(locz)/res)+1;

        % do the split method if needed
        % this portion also does not apply for plasmids
        if dosplit && havebio
            % biomass portion
            agrid(ispec,1).data(i,j,k) = agrid(ispec,1).data(i,j,k) + drow(locb)/vol;
        
            % only add inert mass if we actually have it
			if havenrt
				% inert biomass portion
				agrid(ispec,2).data(i,j,k) = agrid(ispec,2).data(i,j,k) + drow(loci)/vol;
			end

            % only add capsule mass if we actually have it
            if havecap
                % capsule portion
                agrid(ispec,3).data(i,j,k) = agrid(ispec,3).data(i,j,k) + drow(locc)/vol;
            end
        else
            % otherwise we keep all the biomass together
            agrid(ispec,1).data(i,j,k) = agrid(ispec,1).data(i,j,k) + drow(locb)/vol;

            % only add inert mass if we actually have it
			if havenrt
				% inert biomass portion
				agrid(ispec,1).data(i,j,k) = agrid(ispec,1).data(i,j,k) + drow(loci)/vol;
			end

            % only add capsule mass if we actually have it
            if havecap
                % capsule portion
                agrid(ispec,1).data(i,j,k) = agrid(ispec,1).data(i,j,k) + drow(locc)/vol;
            end
        end
    end
end

