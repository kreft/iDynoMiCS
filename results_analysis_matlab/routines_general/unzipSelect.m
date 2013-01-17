function test=unzipSelect(zipFilename,targetFile)

import java.io.*;
import java.util.zip.ZipFile;
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

% This InterruptibleStreamCopier is unsupported and may change without notice.
interruptibleStreamCopier = ...
    InterruptibleStreamCopier.getInterruptibleStreamCopier;


zipFile = ZipFile(zipFilename);
outputDirectory = pwd;

% Inflate all entries.
enumeration = zipFile.entries;
test=false;

while (enumeration.hasMoreElements & ~test)
    zipEntry = enumeration.nextElement;
    if strcmp(zipEntry.getName , targetFile)
        outputName = fullfile(outputDirectory,char(targetFile));
        test=true;
    end
end

if test
    file = java.io.File(outputName);
    parent = File(file.getParent);
    parent.mkdirs;
    try
        fileOutputStream = java.io.FileOutputStream(file);
    catch
        zipFile.close;
        file.close;
        error('Could not create "%s".',outputName);
    end
    % Extract entry to output stream.
    inputStream = zipFile.getInputStream(zipEntry);
    interruptibleStreamCopier.copyStream(inputStream,fileOutputStream);
    % Close streams.
    fileOutputStream.close;
    inputStream.close;
else
    %disp('File not found');
    fprintf('File "%s" not found.\n',targetFile);
end


% Close zip.
zipFile.close;
