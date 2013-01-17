function displaySolute(s)

res = s(1).resolution;
is3D = length(size(s(1).data))==3;

%read solute data and display

th=mean(mean(s(end).data,3),2);
thMax = min(find(th==0));

j = 0;
h1=figure;

for i=1:length(s)
    D   = s(i).data;

    if(length(unique(s(i).data))<=1)
        continue;
    end

    if(j==4)
        h1=figure;
        j=0;
    end

    if is3D
        D   = smooth3(D,'box',3);
    end
    %Plot profile
    subplot(4,3,3*j+1);
    dispProfile(D,s(i).resolution,is3D);
    ylabel(s(i).name);
    
    subplot(4,3,3*j+2);
    dispSurface(D,s(i).resolution,is3D,thMax);

    subplot(4,3,3*j+3);
    dispFlatProfile(D,s(i).resolution,is3D);

    j = j+1;
end

%______________________________________________________________________
function dispProfile(D,res,is3D)

if is3D
	%D = shiftdim(D,1);
    D   = mean(D,3);
end
x = ((1:size(D,2)) - .5)*res;
y = ((1:size(D,1)) - .5)*res;
[X,Y] = meshgrid(x,y);

contourf(X,Y,D);
colorbar


%______________________________________________________________________
function dispSurface(D,res,is3D,thMax)
if is3D
    D   = squeeze(mean(D(1:thMax,:,:),1));
end
x = ((1:size(D,2)) - .5)*res;
y = ((1:size(D,1)) - .5)*res;
[X,Y] = meshgrid(x,y);

contourf(X,Y,D);
colorbar

%______________________________________________________________________
function dispFlatProfile(D,res,is3D)
if is3D
    D   = squeeze(mean(D,3));
end

profile = mean(D,2);
xAxe = [0:length(profile)-1]*res;
plot(xAxe,profile);

%subplot(1,2,2);
%[X,Z] = meshgrid(z,y);
%contourf(X,Z,squeeze(mean(D,2)))
%[X,Y,Z] = meshgrid(x,y,z);
%slice(X,Y,Z,D,[],[10 20 30 40],[80])
%colorbar
