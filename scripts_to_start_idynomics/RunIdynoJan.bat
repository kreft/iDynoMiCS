@echo Off
rem Put this file into the directory with the protocol file to start iDynoMiCS
rem Class files can be in a bin directory or src directory, depending on your settings of Eclipse
rem Jan Kreft j.kreft@bham.ac.uk

rem Set the default xmin and xmax values, and the absolute path to the idynomics folder
SET XMIN=1000m
SET XMAX=1000m
SET NPATH=D:\data\models\eclipse_workspace\iDynoMiCS_June2012
Echo %NPATH%

rem Set the absolute classpath
SET CLASSPATH=%NPATH%\bin
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jcommon-1.0.12.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jfreechart-1.0.9.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\jdom.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\truezip-6.jar
SET CLASSPATH=%CLASSPATH%;%NPATH%\src\lib\Jama-1.0.2.jar

rem Call the program
java -Xms%XMIN% -Xmx%XMAX% -cp %CLASSPATH% idyno.Idynomics %1