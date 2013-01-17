function plotTime(type)
% plotTime(type)
%
% this plots different data over time
%
% for performance reasons, it's best to manually unzip the agent_Sum.zip and
% env_Sum.zip files before running this command
%
% possible arguments:
%
% 'bulks': plots the bulk concentrations over time
% 'bulk_rates': plots the rate of change of the bulks over time
% 'agars': plots the mean agar concentrations over time
% 'N': plots the number of each particle over time
% 'biomass': plots biomass of each species over time
% 'biomass_split': plots biomass of each split species over time
% 'interface': plots the thickness over time
% 'simpsonIndex': plots the diversity index
% 'simpsonIndex_noeps': plots the diversity index but ignores EPS particles
%

if nargin < 1
	help plotTime;
	return;
end

fs = 12;
lw = 1;
fs = 24;
lw = 3;

if strcmp(type,'bulks') || strcmp(type,'bulk_rates')
    iters = getListOfIterates('env_Sum');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping env_Sum.zip\n');
		unzipfiles('env_Sum.zip');
		iters = getListOfIterates('env_Sum');
	end
    
    % first get setup information
    [b,t] = getBulkData(iters(1));
    nBulks = b(1,1).number_bulks;
    nSolutes = b(1,1).number_solutes;

    bulkNames = cell(1,nBulks);
    for i=1:nBulks
        bulkNames{i} = strrep(b(i,1).bulk_name,'_','\_');
    end

    soluteNames = cell(1,nSolutes);
    for i=1:nSolutes
        soluteNames{i} = strrep(b(1,i).solute_name,'_','\_');
    end
    
    % this is where the data will be stored
    sVals = cell(nBulks,1);
    if strcmp(type,'bulks')
        sUnit = b(1).solute_unit;
    else
        sUnit = b(1).uptake_unit;
    end

    % now get the data
    for ib=1:nBulks
        sVals{ib} = zeros(length(iters),nSolutes);
    end
    tVals = zeros(length(iters),1);
    
    for n=1:length(iters)
        [b,t] = getBulkData(iters(n));
        for ib=1:nBulks
            for is=1:nSolutes
                if strcmp(type,'bulks')
                    sVals{ib}(n,is) = b(ib,is).solute_data;
                else
                    sVals{ib}(n,is) = b(ib,is).uptake_data;
                end
            end
        end
        tVals(n) = t;
    end

    % output results
    for ib=1:nBulks
        figure;
        plot(tVals,sVals{ib},'LineWidth',lw);
        lhand = legend(soluteNames);
        set(lhand,'FontSize',fs);
        xlabel('Time [hrs]','FontSize',fs);

        if strcmp(type,'bulks')
            ylabel(sprintf('Concentration [%s]',sUnit),'FontSize',fs);
        else
            ylabel(sprintf('Concentration Rate [%s]',sUnit),'FontSize',fs);
        end
        title(bulkNames{ib},'FontSize',fs);
        set(gca,'FontSize',fs);

		fprintf('Final values:\n');
		disp(sVals{ib}(end,:));
    end

elseif strcmp(type,'agars')
    iters = getListOfIterates('env_Sum');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping env_Sum.zip\n');
		unzipfiles('env_Sum.zip');
		iters = getListOfIterates('env_Sum');
	end
    
    % first get setup information
    [a,t] = getAgarData(iters(1));
    nAgars = a(1,1).number_agars;
    nSolutes = a(1,1).number_solutes;

    agarNames = cell(1,nAgars);
    for i=1:nAgars
        agarNames{i} = strrep(a(i,1).agar_name,'_','\_');
    end

    soluteNames = cell(1,nSolutes);
    for i=1:nSolutes
        soluteNames{i} = strrep(a(1,i).solute_name,'_','\_');
    end
    
    % this is where the data will be stored
    sVals = cell(nAgars,1);
	sUnit = a(1).solute_unit;

    % now get the data
    for ia=1:nAgars
        sVals{ia} = zeros(length(iters),nSolutes);
    end
    tVals = zeros(length(iters),1);
    
    for n=1:length(iters)
        [a,t] = getAgarData(iters(n));
        for ia=1:nAgars
            for is=1:nSolutes
				sVals{ia}(n,is) = mean(mean(a(ia,is).solute_data));
            end
        end
        tVals(n) = t;
    end

    % output results
    for ia=1:nAgars
        figure;
        plot(tVals,sVals{ia},'LineWidth',lw);
        lhand = legend(soluteNames);
        set(lhand,'FontSize',fs);
        xlabel('Time [hrs]','FontSize',fs);
		ylabel(sprintf('Concentration [%s]',sUnit),'FontSize',fs);
        title(agarNames{ia},'FontSize',fs);
        set(gca,'FontSize',fs);

		fprintf('Final values:\n');
		disp(sVals{ia}(end,:));
    end

elseif strcmp(type,'biomass') || strcmp(type,'biomass_split')
	% get a sum of the biomass amount over time
    iters = getListOfIterates('agent_State');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping agent_State.zip\n');
		unzipfiles('agent_State.zip');
		iters = getListOfIterates('agent_State');
	end

    tVals = zeros(length(iters),1);
    BVals = repmat(tVals,1,20);
	maxloc = 1;

    for n=1:length(iters)
		% read in the agent data to a grid, and as split if desired
		if strcmp(type,'biomass')
			[a,t] = getAgentData(iters(n),'Grid');
		else
			[a,t] = getAgentData(iters(n),'SplitGrid');
		end

		nspec = size(a,1);
		ntype = size(a,2);
		for i=1:nspec
			for j=1:ntype
				loc = (i-1)*ntype + j;
				BVals(n,loc) = sum(sum(calcBiofilmData(a(i,j),'biomass')));
				maxloc = max(loc,maxloc);
			end
		end
        tVals(n) = t;
    end
	BVals = BVals(:,1:maxloc);

	% finally get agent names sorted correctly
	aname = {};
	for i=1:nspec
		for j=1:ntype
			loc = (i-1)*ntype + j;
			aname{loc} = a(i,j).name;
		end
	end

	figure;
    plot(tVals,BVals,'LineWidth',lw);
    xlabel('Time [hrs]','FontSize',fs);
    ylabel('Biomass Amount [g/L]','FontSize',fs);
	lhand = legend(aname);
    set(lhand,'FontSize',fs);
    set(gca,'FontSize',fs);

elseif strcmp(type,'simpsonIndex') || strcmp(type,'simpsonIndex_noeps')
    % TODO: right now this treats inert particles as equal to active
    % particles, but they should probably be excluded in finding D
    
    iters = getListOfIterates('agent_Sum');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping agent_Sum.zip\n');
		unzipfiles('agent_Sum.zip');
		iters = getListOfIterates('agent_Sum');
	end
    
    tVals = zeros(length(iters),1);
    DVals = tVals;

    for n=1:length(iters)
		[a,t] = getAgentData(iters(n),'Sum');
        DVals(n) = calcAgentData(a,type);
        tVals(n) = t;
    end

	figure;
    plot(tVals,DVals,'LineWidth',lw);
    xlabel('Time [hrs]','FontSize',fs);
    ylabel('Diversity Index','FontSize',fs);
    set(gca,'FontSize',fs);

elseif strcmp(type,'N')
    iters = getListOfIterates('agent_Sum');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping agent_Sum.zip\n');
		unzipfiles('agent_Sum.zip');
		iters = getListOfIterates('agent_Sum');
	end

    % first get labels
    [a,t] = getAgentData(iters(1),'Sum');
    nSpecies = numel(a);
    labels = cell(nSpecies,1);
    for i=1:nSpecies
        labels{i} = strrep(a(i).name,'_','\_');
    end

    % get the data
    tVals = zeros(length(iters),1);
    nVals = zeros(length(iters),nSpecies);
    for n=1:length(iters)
        [a,t] = getAgentData(iters(n),'Sum');
        for i=1:nSpecies
            nVals(n,i) = a(i).data(1);
        end
        tVals(n) = t;
    end

    % output results
	figure;
    plot(tVals,nVals,'LineWidth',lw);
    lhand = legend(labels);
    set(lhand,'FontSize',fs);
    xlabel('Time [hrs]','FontSize',fs);
    ylabel('Number of Particles','FontSize',fs);
    set(gca,'FontSize',fs);

elseif strcmp(type,'interface')
    % get the interface data from the first species
    
    iters = getListOfIterates('env_Sum');
	if numel(iters)==0
		fprintf('There are no iterates to plot; unzipping env_Sum.zip\n');
		unzipfiles('env_Sum.zip');
		iters = getListOfIterates('env_Sum');
	end
    
    tVals = zeros(length(iters),1);
    iVals = zeros(length(iters),2);

    for n=1:length(iters)
        [d,t,h] = getBulkData(iters(n));
		if h(1) > 0
			% h has height data as (mean, stddev, max)
			iVals(n,1) = h(3);
			iVals(n,2) = h(1);
		else
			% old method that requires looking at biomass distributions to infer biofilm height
			fprintf('Using old method to get location for iterate %i, time %g.\n',n,t);
			[s,a,r,t,h] = getEnvData(n);
			idata = calcBiofilmData(a,'interface');
			iVals(n,1) = max(idata);
			iVals(n,2) = mean(idata);
		end
        tVals(n) = t;
    end

    % output results
	figure;
    plot(tVals,iVals,'LineWidth',lw);
    lhand = legend('Maximum','Average','Location','NorthWest');
    set(lhand,'FontSize',fs);
    xlabel('Time [hrs]','FontSize',fs);
    ylabel('Biofilm Thickness [\mum]','FontSize',fs);
    set(gca,'FontSize',fs);
end


