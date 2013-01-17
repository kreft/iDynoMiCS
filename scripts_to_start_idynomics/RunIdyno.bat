@echo Off
rem Multi_run_idyno.bat v1.0
rem Based on run_idyno.bat.old
rem Author: Edd Miles - edd@edd-miles.com
rem Script to start the idynomics software from the command line under Windows.
rem Three (optional) command switches may be provided:
rem *DEPRECIATED* -i for the protocol files to be worked on *DEPRECIATED*
rem Input files no longer need to be preceeded by -i
rem -min for the java -Xms
rem -max from the java -Xmx variable
rem For ease of use, it is advised to add the idynomics folder to your PATH variable

rem Set the default xmin and xmax values, and the absolute path to the idynomics folder
SET XMIN=1000m
SET XMAX=1000m
SET NPATH=D:\data\models\eclipse_workspace\iDynoMiCS_June2012
Echo %NPATH%

IF EXIST Tfile.txt del Tfile.txt
rem Hacky method of allowing command line switches in bat files. If anyone knows
rem a more elegant solution, please email me at edd@edd-miles.com
:LoopStart
IF "%1" == "-min" SET XMIN=%2 & SHIFT & GOTO LoopBottom
IF "%1" == "-max" SET XMAX=%2 & SHIFT & GOTO LoopBottom
IF "%1" == "-i" SET SET iFILE=%2 & SHIFT & ( Echo Warning: Use of -i is depreciated ) & GOTO FileAppend
IF "%1" == "" GOTO LoopEnd
SET iFILE=%1 & GOTO FileAppend

:FileAppend
IF EXIST %iFILE% GOTO MidFileAppend
Echo "NO SUCH FILE %iFILE%"
exit
:MidFileAppend
Echo %iFILE% >> Tfile.txt || ( echo "FAILED TO WRITE TEMP FILE" & exit )

:LoopBottom
SHIFT
GOTO LoopStart
:LoopEnd

rem Set the absolute classpath
SET CLASSPATH=%NPATH%\bin
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jcommon-1.0.12.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jfreechart-1.0.9.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jdom.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\truezip-6.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\Jama-1.0.2.jar


rem call the program
for /f "tokens=* delims= " %%a in (Tfile.txt) do ( java -Xms%XMIN% -Xmx%XMAX% -cp %CLASSPATH% idyno.Idynomics %%a )
DEL Tfile.txt
