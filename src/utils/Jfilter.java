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

/**
 * \brief Used in simulation initialisation to ensure onky XML files are read as protocol files and previous results files are excluded
 * 
 * Used in simulation initialisation to ensure onky XML files are read as protocol files and previous results files are excluded. This 
 * can occur when a simulation is being restarted
 */
public class Jfilter implements FilenameFilter 
{

	/**
	 * File extension by which filter will be applied
	 */
	String _extension;

	/**
	 * \brief Initialised the file extension filter and sets the extension to the specified argument
	 * 
	 * Initialised the file extension filter and sets the extension to the specified argument
	 * 
	 * @param extension	The file extension that should be filtered
	 */
	public Jfilter(String extension) 
	{
		_extension = "."+extension;
	}

	/**
	 * \brief Used to determine if a file should be excluded from the filter or whether the file name is acceptable
	 * 
	 * Used to determine if a file should be excluded from the filter or whether the file name is acceptable
	 * 
	 * @param directory	The directory of files being filtered
	 * @param filename	The filename being checked against the filter
	 * @return Boolean value noting whether or not the file ends with the extension
	 */
	public boolean accept(File directory, String filename) 
	{

		if (filename.endsWith(_extension)) return true;
		return false;

	}
}
