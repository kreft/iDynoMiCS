function plotProfile(nFile,yloc,xmax)
% plotProfile(nFile,yloc,xmax)
%
% this plots the profile of biomass and solutes as a function of
% x (height) for a given horizontal position (yloc)
% 
% xmax is optional, but allows you to limit the vertical space included
% set xmax to 'edge' to limit the plot to the within-biofilm region
%
% alternatively, set yloc to 'avg' or 'std' to plot the average
% or standard deviation of values across the domain width
%
% NOTE: the plot labels assume all components have the same units
%

if nargin < 1
	help plotProfile;
	return;
end

if nargin < 2
    yloc = 'avg';
end

if nargin < 3
    xmax = 10000;
end

if ~strcmp(xmax,'edge') && ~isnumeric(xmax)
    fprintf('xmax argument is invalid.\n');
    xmax = 10000;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we don't need to read in the agent data directly
% because we will use the averages given in the env_State file
%[a,time] = getAgentData(nFile,'State');
%nspec = numel(a);

[s,a,r,time] = getEnvData(nFile);

nsols = numel(s);
nspec = numel(a);
nreac = numel(r);

% now replace the agent data with that read in from the agent file, which
% allows us to split apart the active biomass, inert biomass, and EPS
% portions
a = getAgentData(nFile,'SplitGrid');
nspec = size(a,1);
ntype = size(a,2);


% make names legible in plots (need to escape the underscore)
for i=1:nsols
    s(i).name = strrep(s(i).name,'_','\_');
end
for i=1:nspec
    for j=1:ntype
        a(i,j).name = strrep(a(i,j).name,'_','\_');
    end
end
for i=1:nreac
    r(i).name = strrep(r(i).name,'_','\_');
end

% at this point we have:
% s structure for solutes; all nsols entries are the solutes
% a structure for agents; all nspec entries are for agents; biomass, EPS,
%      and inerts are separated
% r structure for reactions; all nreac entries are for reactions

fprintf('Plotting results for time %g hours.\n',time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spatial information
res = s(1).resolution;
nx = s(1).nI;
ny =s(1).nJ;
x = (((1:nx) - .5)*res);
y = (((1:ny) - .5)*res);

% spatial information for agents (uses the agentGrid rather than solute grid)
ares = a(1,1).resolution;
anx = a(1).nI;
any =a(1).nJ;
ax = (((1:anx) - .5)*ares);
ay = (((1:any) - .5)*ares);

% make sure yloc is not outside the bounds
if isnumeric(yloc)
    yloc = max(y(1),yloc);
    yloc = min(y(ny),yloc);
    fprintf('Plotting profile at location %i um.\n',yloc);
end

%[X,Y] = meshgrid(x,y);
%contourf(X,Y,(s(1).data)');

sdata = zeros(nsols,nx);
adata = zeros(nspec*ntype,anx);
rdata = zeros(nreac,nx);

if strcmp(yloc,'avg') || strcmp(yloc,'std')
    % calculate mean value
    for i=1:nsols
        sdata(i,:) = mean(s(i).data,2);
    end
    for i=1:nspec
        for j=1:ntype
            loc = (i-1)*ntype + j;
            adata(loc,:) = mean(a(i,j).data,2);
            
        end
    end
    for i=1:nreac
        rdata(i,:) = mean(r(i).data,2);
    end
end
if strcmp(yloc,'std')
    % also calculate standard deviation if needed
    for i=1:nsols
        sstdv(i,:) = std(s(i).data',1);
    end
    for i=1:nspec
        for j=1:ntype
            loc = (i-1)*ntype + j;
            astdv(loc,:) = std(a(i,j).data',1);
        end
    end
    for i=1:nreac
        rstdv(i,:) = std(r(i).data',1);
    end
end
if isnumeric(yloc)
    % calculate interpolated value when location is given
    for i=1:nsols
        sdata(i,:) = interp2(x,y,s(i).data',x,yloc);
    end
    for i=1:nspec
        for j=1:ntype
            loc = (i-1)*ntype + j;
            adata(loc,:) = interp2(ax,ay,a(i,j).data',ax,yloc);
        end
    end
    for i=1:nreac
        rdata(i,:) = interp2(x,y,r(i).data',x,yloc);
    end
end

% location of the interface
if isnumeric(yloc)
	% calculate interpolated value
    xinter = interp1(ay,calcBiofilmData(a,'interface'),yloc);
else
	% assume we need the mean value if location isn't specified
	%xinter = mean(calcBiofilmData(a,'interface'));
	[bdata,btime,heights] = getBulkData(nFile);
	xinter = heights(1);
end

% finally get agent names sorted correctly
aname = {};
for i=1:nspec
	for j=1:ntype
		loc = (i-1)*ntype + j;
		aname{loc} = a(i,j).name;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 3;
fs = 24;

maintitle = sprintf('Profile at %g hours',time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solutes
figure;
if strcmp(yloc,'std')
    xs = repmat(x,nsols,1);
    errorbar(xs',sdata',sstdv','LineWidth',lw);
else
    plot(x,sdata,'LineWidth',lw);
end
ylimsave = ylim;
hold on;
plot([xinter,xinter],[-10000,10000],'--r','LineWidth',lw);

if strcmp(xmax,'edge')
    xmax = xinter;
else
    xmax = min(max(x),xmax);
end
xlim([min(x) xmax]);
ylim(ylimsave);

title(maintitle,'FontSize',fs);
xlabel('X [\mum]','FontSize',fs);
ylabel(sprintf('Solute Concentration [%s]',s(1).unit),'FontSize',fs);
lhand = legend(s.name);
set(lhand,'FontSize',fs);
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biomass

% (this little bit allows comparison with previous routine by just totaling
%  the amount of biomass for each species, like we had previously)
% newa = zeros(nspec,anx);
% for i=1:nspec
%     for j=1:ntype
%         newa(i,:) = newa(i,:) + adata((i-1)*ntype+j,:);
%     end
% end
% adata = newa;

figure;
if strcmp(yloc,'std')
    xs = repmat(ax,nspec*ntype,1);
    errorbar(xs',adata',astdv','LineWidth',lw);
else
    plot(ax,adata,'LineWidth',lw);
end

ylimsave = ylim;
hold on;
plot([xinter,xinter],[-10000,10000],'--r','LineWidth',lw);

xlim([min(ax) xmax]);
ylim(ylimsave);

title(maintitle,'FontSize',fs);
xlabel('X [\mum]','FontSize',fs);
ylabel(sprintf('Biomass Concentration [%s]',a(1).unit),'FontSize',fs);
lhand = legend(aname);
set(lhand,'FontSize',fs);
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactions
figure;
if strcmp(yloc,'std')
    xs = repmat(x,nreac,1);
    errorbar(xs',rdata',rstdv','LineWidth',lw);
else
    plot(x,rdata,'LineWidth',lw);
end

ylimsave = ylim;
hold on;
plot([xinter,xinter],[-10000,10000],'--r','LineWidth',lw);

xlim([min(x) xmax]);
ylim(ylimsave);

title(maintitle,'FontSize',fs);
xlabel('X [\mum]','FontSize',fs);
ylabel(sprintf('Reaction Rate [%s]',r(1).unit),'FontSize',fs);
lhand = legend(r.name);
set(lhand,'FontSize',fs);
set(gca,'FontSize',fs);


