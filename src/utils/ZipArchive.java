/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import java.io.*;
import java.util.zip.*;

/**
 * \brief Class to create and update zip archives that contain simulation result files
 * 
 * iDynoMiCS simulations can produce a large number of result files, depending on the output period specified. With each timestep these 
 * files are added to a zip archive to keep these files together. This class contains methods to create and update these files
 *
 */
public class ZipArchive 
{

	/**
	 * Memory buffer for reading result file information
	 */
	static final int        BUFFER = 2048;
	
	/**
	 * Output stream used to create the zip file of results
	 */
	private ZipOutputStream out;
	
	/**
	 * The name assigned to the zip file
	 */
	private String          _zipFileName;

	/**
	 * \brief Constructor to open a zip file output stream and set the zip file name
	 * 
	 * Constructor to open a zip file output stream and set the zip file name
	 * 
	 * @param zipFileName	The name assigned to this zip file of simulation results
	 * @throws IOException	Exception thrown if this stream cannot be opened
	 */
	public ZipArchive(String zipFileName) throws IOException 
	{
		// create file
		_zipFileName = zipFileName;
		out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(_zipFileName)));
	}

	/**
	 * \brief Adds a simulation results file to its respective zip archive and then deletes the result file
	 * 
	 * Add a simulation results file to its respective zip archive. The file itself is deleted after this has occurred
	 * 
	 * @param f	The results file to be added to the zip archive
	 * @throws IOException	Exception thrown if this zip file cannot be opened
	 */
	public final void addToZipArchiveAndDelete(File f) throws IOException 
	{
		int count;
		byte data[] = new byte[BUFFER];

		BufferedInputStream origin = new BufferedInputStream(new FileInputStream(f), BUFFER);
		out.putNextEntry(new ZipEntry(f.getName()));

		while ((count = origin.read(data, 0, BUFFER))!=-1) 
		{
			out.write(data, 0, count);
		}
		out.closeEntry();
		origin.close();
		f.delete();
	}

	/**
	 * \brief Adds a simulation results file to its respective zip archive, and keeps the result file under a specified name
	 * 
	 * Add a simulation results file to its respective zip archive. In this version of the method the result file is retained
	 * 
	 * @param f	The results file to be added to the zip archive
	 * @param fNew	The new name to be assigned to this results file
	 * @throws IOException	Exception thrown if this stream cannot be opened
	 */
	public final void addToZipArchiveAndKeepLast(File f, File fNew) throws IOException 
	{
		int count;
		byte data[] = new byte[BUFFER];

		BufferedInputStream origin = new BufferedInputStream(new FileInputStream(f), BUFFER);
		out.putNextEntry(new ZipEntry(f.getName()));

		while ((count = origin.read(data, 0, BUFFER))!=-1) {
			out.write(data, 0, count);
		}
		origin.close();
		f.renameTo(fNew);
	}

	/**
	 * \brief	Closes the zip file results archive
	 * 
	 * Closes the zip file results archive
	 * @throws IOException	Exception thrown if the file cannot be closed
	 */
	public void close() throws IOException {
		out.close();
	}
}
