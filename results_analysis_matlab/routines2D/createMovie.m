function createMovie(movieName,useExisting,docrop)
% createMovie(movieName,useExisting,docrop)
%
% This creates a movie by creating the png files from
% the it(*).pov files.
%
% Set 'useExisting' to 1 to skip creating the image files and
% instead use what was previously rendered
%
% The docrop argument is optional. If 'docrop' is 0 then
% the images will NOT be cropped. Default value is 1.
%
%


if nargin < 2
    useExisting = 0;
end
if nargin < 3
    docrop = 1;
end

% first create the needed image files
if ~useExisting
    createPicsNew('all',docrop);
end

% now make the movie based on the type of output file

outtype = movieName(end-2:end);


if strcmp(outtype,'avi')
	% create the movie using Matlab routines

	iters = getListOfIterates('it','png');
	mov = avifile(movieName);
	for it=1:length(iters)
		theim = imread(sprintf('it(%04i).png',iters(it)),'BackgroundColor',[1,1,1]);
		mov = addframe(mov,im2frame(theim));
	end
	mov = close(mov);

elseif strcmp(outtype,'gif')
	% create a movie using ImageMagick

	imPATH = getProgramPathNew('ImageMagick');
    dos(sprintf('"%s\\convert.exe" -delay 10 *.png %s',imPATH,movieName));
  
    
elseif strcmp(outtype,'mpg')
	% create a movie using ImageMagick
	imPATH = getProgramPathNew('ImageMagick');
	dos(sprintf('"%s\\convert.exe" -delay 10 *.png %s',imPATH,movieName));

	%% create a movie using ffmpeg
	%ffmpegPATH = getProgramPath('ffmpeg');
	%dos(sprintf('%s\\ffmpeg -v -l -i it(%s).png %s',ffmpegPATH,'%04d',movieName));
end


