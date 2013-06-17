/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import idyno.SimTimer;
import java.io.BufferedOutputStream;
import de.schlichtherle.io.File;
import de.schlichtherle.io.FileOutputStream;

/**
 * \brief Class used to create auto-zipping result files. At each iteration, the file is added to an archive.
 * 
 * Class used to create auto-zipping result files. At each iteration, the file is added to an archive. The last iteration is still 
 * available with the syntax "resultFileName(last)"
 *   
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */
public class ResultFile 
{
	/**
	 * Result file archive (zip file) in which this result file will be stored
	 */
	private File	_vArchive;
	
	/**
	 * Name of the new results file being created at the required output period
	 */
	private File	newFile;
	
	/**
	 * Name of the archive file (.zip) that the new results file will be included within
	 */
	private File	archiveFile;
	
	/**
	 * Simulation iteration at which this file is being written
	 */
	private int                  _fileIndex;
	
	/**
	 * Output stream used to write a results file
	 */
	private FileOutputStream     _out;

	/**
	 * Name assigned to this result file
	 */
	private String               _prefix;
	
	/**
	 * String noting the directory this result file should be written to
	 */
	private String               _dir;

	/**
	 * Buffered output stram used to write simulation statistics to file
	 */
	private BufferedOutputStream buffer;
	
	/**
	 * Buffer to hold the information / statistics that are being written to the results file
	 */
	private StringBuffer         value;

	/**
	 * \brief Creates a result file within a specified folder, of a particular name, and at a set simulation iteration
	 * 
	 * Creates a result file within a specified folder, of a particular name, and at a set simulation iteration
	 * 
	 * @param outPath	directory where result files will be saved
	 * @param fileName	name of the resultFile
	 * @param iter	first iterate to start writing to (added bvm 26.1.2009)
	 */
	public ResultFile(String outPath, String fileName, int iter) 
	{
		// Set directory and filename
		_prefix = fileName;
		_dir = outPath+File.separator;
		_fileIndex = iter;

		// Create the archive file
		_vArchive = new File(_dir+_prefix+".zip");
		if (!_vArchive.exists()) _vArchive.mkdir();		
	}

	/**
	 * \brief Creates a result file for the current iteration
	 * 
	 * Create a resultFile for the current iteration
	 * 
	 * @throws Exception	Generated exception should read errors occur
	 * 
	 * @param iter	The current simulation iteration (added by BVM Jan 2009)
	 */
	public void openFile(int iter) 
	{
		try 
		{
			// bvm added 26.1.2009: use simulation iterate for file name 
			_fileIndex = iter;

			newFile = new File(_dir+"lastIter"+File.separator+_prefix+"(last).xml");
			archiveFile = new File(_dir+_prefix+".zip"+File.separator
									+_prefix+"("+_fileIndex+").xml");

			// Create the streams to write in the file
			_out = new FileOutputStream(newFile);
			buffer = new BufferedOutputStream(_out);

			// Build the main markup
			// bvm 26.1.2009: added output of iterate as well as time
			value = new StringBuffer("<idynomics>\n <simulation iterate=\"");
			value.append(SimTimer.getCurrentIter());
			value.append("\" time=\"");
			value.append(SimTimer.getCurrentTime());
			value.append("\" unit=\"hour\">\n");
			buffer.write(value.toString().getBytes());
			
		} catch (Exception e) 
		{
			LogFile.writeError("Unable to open result file", "ResultFile.openFile()");
		}
	}

	/**
	 * \brief Add text to an existing resultFile
	 * 
	 * Add text to an existing resultFile
	 * 	 * 
	 * @param text	The text that should be appended to this results file
	 */
	public void write(String text) {
		try {
			buffer.write(text.getBytes());
		} catch (Exception e) {
			LogFile.writeError("Unable to write", "ResultFile.write()");
		}
	}

	/**
	 * \brief Closes the resultFile and adds the closing mark-up to have a well-formed XML file
	 * 
	 * Closes the resultFile and adds the closing mark-up to have a well-formed XML file
	 */
	public void closeFile() {
		try {
			// Close the markup
			buffer.write("\n</simulation>\n</idynomics>".getBytes());

			// Close the file
			buffer.close();
			_out.close();

			// Add the resultFile to the archive
			newFile.copyTo(archiveFile);
			File.update(_vArchive);

		} catch (Exception e) {
			LogFile.writeError("Unable to close or archive", "ResultFile.closeFile()");
		}
	}

	/**
	 * \brief Static function invoked to copy a File of a given name to a specified target
	 * 
	 * Static function invoked to copy a File of a given name to a specified target
	 * 
	 * @param sourceName	The name of the results file to copy
	 * @param targetName	The name that should be given to the copy
	 * @return Boolean noting the result of the copy procedure
	 */
	public static boolean copyFile(String sourceName, String targetName) {		
		File source = new File(sourceName);
		File target = new File(targetName);
		return source.copyTo(target);
	}
}
