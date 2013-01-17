function makePic(nBeg,nEnd,nStep)

POVPATH='c:\apps\science\POVRay\bin\';
QuietPATH='c:\apps\science\QuietPOV\';

%Start POV-Ray
eval(['!',POVPATH,'pvengine &']);
pause(1);

for iFile=nBeg:nStep:nEnd
    fName = ['it(',num2str(iFile),').pov']
    unzipSelect('povray.zip',fName);
    unzipSelect('povray.zip','sceneheader.inc');
    unzipSelect('povray.zip','scenefooter.inc');
    %eval(['!',QuietPATH,'quietpov -start ',fName]);
    dos([QuietPATH,'quietpov -start ',fName]);
    %delete(fName);
end
