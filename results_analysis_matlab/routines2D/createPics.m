function createPics(nFile,docrop)
% createPics(nFile,docrop)
%
% This function will create png images from pov files.
% nFile may be set to a file number, 'last', or 'all'.
%
% The docrop argument is optional. If 'docrop' is 1 then
% the image will be cropped to within a small border
% around the content.
%

if nargin < 1
	help createPics;
	return;
end
if nargin < 2
    docrop = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get paths to programs needed here
global POVPATH QuietPATH imPATH;

POVPATH   = getProgramPath('POV-Ray');
QuietPATH = getProgramPath('QuietPOV');
imPATH    = getProgramPath('ImageMagick');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start POV-Ray
%eval(['!',POVPATH,'pvengine &']);
dos(sprintf('"%s\\pvengine.exe" &',POVPATH));

pause(1);

% get header and footer files used everywhere
if ~exist('sceneheader.inc','file')
    %unzipSelect('povray.zip','sceneheader.inc');
	copyfile('lastIter/sceneheader.inc','.');
end
if ~exist('scenefooter.inc','file')
    %unzipSelect('povray.zip','scenefooter.inc');
	copyfile('lastIter/scenefooter.inc','.');
end

% now the treatment differs based on input arguments
if isnumeric(nFile)
	% if a particular iterate was specified
    theFile = sprintf('it(%i).pov',nFile);
    if ~exist(theFile,'file')
        unzipSelect('povray.zip',theFile);
    end
    theImage = renderPOV(nFile,theFile);
	if docrop
		cropImages(theImage);
	end
end
if strcmp(nFile,'last')
	% if we want to render the last iterate
    theFile = 'lastIter/it(last).pov';
    theImage = renderPOV(nFile,theFile);
	if docrop
		cropImages(theImage);
	end
end
if strcmp(nFile,'all')
	
    % if we want to render all iterates

    if exist('povray.zip','file')
    fprintf('Unzipping povray files... \n');
    unzipfiles('povray.zip');
    fprintf('Done.\n');
    end
    
    iters = getListOfIterates('it','pov');
	allImages = {};
    for it=1:length(iters)
        nFile = iters(it);
        theFile = sprintf('it(%i).pov',nFile);
		if ~isnumeric(nFile)
			theImage = 'it(last).png';
		else
			theImage = sprintf('it(%04i).png',nFile);
		end

		if exist(theImage,'file')
			fprintf('Image %s exists; skipping rendering.\n',theImage);
		else
			theImage = renderPOV(nFile,theFile);
		end
		allImages{end+1} = theImage;
    end
	if docrop
		cropImages(allImages);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imgfile = renderPOV(nFile,theFile)

global QuietPATH POVPATH;

fprintf('Rendering file %s\n',theFile);

% do the render and convert the file format
%try using QuietPov, if it doesn't work, skip it and call pvengine directly
%dos([QuietPATH,'quietpov -start ',theFile]);
%dos(sprintf('"%s\\quietpov" -start %s',QuietPATH,theFile));


 dos(sprintf('"%s\\pvengine.exe" /RENDER %s', POVPATH, theFile));
 pause(5)
    
bmp2png(theFile);

if isnumeric(nFile)
    % now pad the file name with zeros so that it sorts correctly
    oldName = sprintf('it(%i).png',nFile);
    newName = sprintf('it(%04i).png',nFile);
	if ~strcmp(oldName,newName)
		movefile(oldName,newName);
	end
else
    newName = 'it(last).png';
	% moves the image file into the base directory
	movefile('lastIter/it(last).png','.');
end

imgfile = newName;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bmp2png(povfile)
% this converts the output .bmp file to .png if it doesn't exist

global imPATH;

bmpfile = strrep(povfile,'pov','bmp');
pngfile = strrep(povfile,'pov','png');

if ~exist(pngfile,'file')
	% now call the convert command to convert the image
	dos(sprintf('"%s\\convert" %s  %s',...
		imPATH,bmpfile,pngfile));
	delete(bmpfile);
end

