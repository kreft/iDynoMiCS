function result = calcAgentData(adata,theval)
% result = calcAgentData(adata,theval)
%
% THESE DO NOT OPERATE SEPARATELY ON SPECIES; IF YOU WANT SEPARATE
% SPECIES DATA THEN ONLY PASS IN THE SINGLE SPECIES
%
% 'number': returns total number of agents in the adata
% 'area': returns area value of all given agents
%         (this is NOT the area taken up by the 'population' b/c
%          there is extra space between agents that this does not include)
% 'simpsonindex': diversity value (0: homogeneous; 1: diverse)
% 'simpsonindex_noeps': diversity value, but ignoring EPS particles
% 'averagesize': returns the average radius of the agents
%
% THE FOLLOWING ASSUME PERIODICITY IN THE HORIZONTAL DIRECTION (Y ONLY)
% 'center': returns center of mass of all given agents in xyz coordinates
% 'stddev': returns std deviation in xyz coordinates
% 'spread': returns maximum spread in xyzcoordinates
%
% (append '_nonper' to the previous routines for non-periodic version: e.g. 'center_nonper')
%
%

nspec = size(adata,1);

result = 0;

theval = lower(theval);

if strcmp(theval,'number')
    result = getNum(adata);

elseif strcmp(theval,'simpsonindex') || strcmp(theval,'simpsonindex_noeps')
	nEach = zeros(nspec,1);

	for i=1:nspec
		if strcmp(theval,'simpsonindex') || ...
				length(findstr('EPS',adata(i).name)) == 0
			% if we're counting all species or it is not an EPS particle
			nEach(i) = getNum(adata(i,1));
		end
	end
	nTot = sum(nEach);

	% now calculate the Simpson Index
	result = 1 - sum(nEach.*(nEach-1))/(nTot*(nTot-1));

elseif strcmp(theval,'area')
    for ispec=1:nspec
        rads = adata(ispec).data(:,adata(ispec).header.radius);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RIGHT NOW THIS ASSUMES A 2D GRID - WILL NEED TO FIX THIS FOR
        % THE GENERAL CASE!!!! TODO TODO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        result = result + pi*sum(rads.*rads);
    end

elseif strcmp(theval,'center') || strcmp(theval,'stddev') || strcmp(theval,'spread')
    [xc,yc,zc] = getCoords(adata);

    % in order to determine these values, use the periodicity of the domain
    % and check offets by integer values of the solute grid to find which
    % offset yields the smallest spread in position values; this will be
    % the offset that clumps the agents the best

	res = adata(1,1).resolution;
	ny = adata(1,1).nJ;

	ioff = 0;

    ymax = res*ny;
    spd = res*ny;
    for i=0:ny-1
        yoff = i*res;
        ycn = yc + (ymax - yoff)*(yc < yoff) - yoff*(yc > yoff);

        ycnm = mean(ycn);
        ycns = std(ycn - ycnm);

        if ycns < spd
            spd = ycns;
            ioff = i;
        end
    end

    % now set up new coordinates based on the offset that yields the lowest
    % spread in values
    yoff = ioff*res;
    yc = yc + (ymax - yoff)*(yc < yoff) - yoff*(yc > yoff);
    % keep the mean within the domain
    ycmean = mod(mean(yc)+yoff,ymax);

    
    % calculate & return the correct value
    if strcmp(theval,'center')
        result = [mean(xc),ycmean,mean(zc)];
    elseif strcmp(theval,'stddev')
        result = [std(xc-mean(xc)),std(yc-ycmean),std(zc-mean(zc))];
    elseif strcmp(theval,'spread')
        % find the spreads based on the middle 90% of agents as a way to
        % ignore outliers (but only if we have enough agents)
        xcs = sort(xc);
        ycs = sort(yc);
        zcs = sort(zc);

        nout = 1;
        nagents = length(xcs);
		if nagents ==0
			result = [0,0,0];
		else
			if nagents > 40
				nout = ceil(0.05*nagents);
			end

			xcmin = xcs(nout);
			xcmax = xcs(nagents+1-nout);
			ycmin = ycs(nout);
			ycmax = ycs(nagents+1-nout);
			zcmin = zcs(nout);
			zcmax = zcs(nagents+1-nout);

			result = [xcmax-xcmin,ycmax-ycmin,zcmax-zcmin];
		end
    end

elseif strcmp(theval,'center_nonper')
    [xc,yc,zc] = getCoords(adata);
    xc = mean(xc);
    yc = mean(yc);
    zc = mean(zc);
    result = [xc,yc,zc];

elseif strcmp(theval,'stddev_nonper')
    [xc,yc,zc] = getCoords(adata);
    xc = std(xc);
    yc = std(yc);
    zc = std(zc);
    result = [xc,yc,zc];

elseif strcmp(theval,'spread_nonper')
    [xc,yc,zc] = getCoords(adata);
    xc = max(xc) - min(xc);
    yc = max(yc) - min(yc);
    zc = max(zc) - min(zc);
    result = [xc,yc,zc];

elseif strcmp(theval,'averagesize')
	% calculate the average size for each species
	result = zeros(nspec,1);
	for ispec=1:nspec
		result(ispec) = mean(adata(ispec).data(:,adata(ispec).header.radius));
	end

else
	% safety catch
	fprintf('\n*** Calculation %s is not known. ***\n',theval);
	fprintf('\t(use help to see options)\n');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getNum(specdata)

n = 0;
for ispec=1:size(specdata,1)
	if isfield(specdata(ispec,1).header,'population')
		n = n + specdata(ispec).data(specdata(ispec,1).header.population);
	else
		n = n + size(specdata(ispec).data,1);
	end
end

%%%%%%%%%%%%%%%%%%%%

function [xc,yc,zc] = getCoords(specdata)

xc = [];
yc = [];
zc = [];

for ispec=1:size(specdata,1)
	% if there are no agents of a species, nothing gets added to the list
	if size(specdata(ispec).data,1) > 0
		xc = [xc; specdata(ispec).data(:,specdata(ispec).header.locationX)];
		yc = [yc; specdata(ispec).data(:,specdata(ispec).header.locationY)];
		zc = [zc; specdata(ispec).data(:,specdata(ispec).header.locationZ)];
	end
end

