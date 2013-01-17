function makeMovie(nD,fName)
ffmpegPATH = 'C:\apps\science\ffmpg\';

%cd povfiles
eval(['!',ffmpegPATH,'ffmpeg -r 2.5 -b 1800 -i it(%0',num2str(nD),'d).png ',fName]);
%cd ..
