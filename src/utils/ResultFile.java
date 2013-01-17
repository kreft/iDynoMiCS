/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *_______________________________________________________
 * ResultFile : class used to create auto-zipping result files
 * At each iteration, the file is added to an archive.
 * The last iteration is still available with the syntax "resultFileName(last)"
 * 
 */

/**
 * @since June 2008
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */

package utils;

import idyno.SimTimer;

import java.io.BufferedOutputStream;
import de.schlichtherle.io.File;
import de.schlichtherle.io.FileOutputStream;

public class ResultFile {

	private File                 _vArchive,newFile, archiveFile;
	private int                  _fileIndex;
	private FileOutputStream     _out;

	private String               _prefix;
	private String               _dir;

	private BufferedOutputStream buffer;
	private StringBuffer         value;

	/* ___________ CONSTRUCTOR _______________________ */

	/**
	 * @param outPath:directory where result files will be saved
	 * @param fileName: name of the resultFile
	 * @param iter: first iterate to start writing to (added bvm 26.1.2009)
	 */
	public ResultFile(String outPath, String fileName, int iter) {
		// Set directory and filename
		_prefix = fileName;
		_dir = outPath+File.separator;
		_fileIndex = iter;

		// Create the archive file
		_vArchive = new File(_dir+_prefix+".zip");
		if (!_vArchive.exists()) _vArchive.mkdir();		
	}

	/**
	 * Create a resultFile for the current iteration
	 * @throws Exception
	 * 
	 * @param iter (bvm added 26.1.2009 to make output more robust)
	 */
	public void openFile(int iter) {
		try {
			// Prepare the File object for the (last) file and the archived one;
			//_fileIndex++;
			
			// bvm added 26.1.2009: use simulation iterate rather than incremental iterate 
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
			
		} catch (Exception e) {
			LogFile.writeError("Unable to open result file", "ResultFile.openFile()");
		}
	}

	/**
	 * Add text to an existing resultFile
	 * @param text
	 */
	public void write(String text) {
		try {
			buffer.write(text.getBytes());
		} catch (Exception e) {
			LogFile.writeError("Unable to write", "ResultFile.write()");
		}
	}

	/**
	 * Close the resultFile (adds the closing mark-up to have a well-formed XML
	 * file
	 * 
	 * @throws Exception
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
	 * Static function invoked to copy a File
	 * @param sourceName
	 * @param targetName
	 * @return the result of the copy procedure
	 */
	public static boolean copyFile(String sourceName, String targetName) {		
		File source = new File(sourceName);
		File target = new File(targetName);
		return source.copyTo(target);
	}
}
