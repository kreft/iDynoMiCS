function thepath = getProgramPath(theprog)
% thepath = getProgramPath(theprog)
% returns the path to a few useful programs
%

theprog = lower(theprog);

if strcmp(theprog,'POV-Ray')
	% install location for POV-Ray
	thepath = 'C:\Program Files (x86)\POV-Ray\v3.6\bin\';

elseif strcmp(theprog,'quietpov')
	% install location for the QuietPOV add-on
	thepath = 'C:\Program Files\POV-Ray for Windows v3.6\guiext\QuietPOV';

elseif strcmp(theprog,'imagemagick')
	% install location for ImageMagick
	thepath = 'C:\Program Files\ImageMagick-6.6.9-Q16';

elseif strcmp(theprog,'ffmpeg')
	% install location for the ffmpeg library
	thepath = 'D:\models\eclipse_workspace\idynomics\software_utilities\ffmpeg';

else
	thepath = '';
end

