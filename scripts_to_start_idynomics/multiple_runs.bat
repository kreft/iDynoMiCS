REM multiple_runs syntax:
REM multiple_runs protocolFile numberOfRuns
REM Example: multiple_runs myProtocol.xml 10

SET CLASSPATH=..\src
SET CLASSPATH=%CLASSPATH%;..\src\lib\jcommon-1.0.12.jar
SET CLASSPATH=%CLASSPATH%;..\src\lib\jfreechart-1.0.9.jar
SET CLASSPATH=%CLASSPATH%;..\src\lib\jdom.jar
SET CLASSPATH=%CLASSPATH%;..\src\lib\truezip-6.jar
SET CLASSPATH=%CLASSPATH%;..\src\lib\Jama-1.0.2.jar

for /l %%n in (1, 1, %2) do java -Xms1000m -Xmx1000m -XX:-UseGCOverheadLimit -cp %CLASSPATH% idyno.Idynomics %1