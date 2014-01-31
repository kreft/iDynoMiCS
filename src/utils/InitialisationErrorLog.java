/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import java.io.File;
import java.io.FileOutputStream;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import idyno.Idynomics;

/**
 * \brief Creates and updates a simulation initialisation error log file for
 * the initialising simulation.
 */
public class InitialisationErrorLog 
{
	/**
	 * Output stream which log messages are written to
	 */
	public static FileOutputStream log;
	
	/**
	 * Format of the date which is used in logging simulation messages
	 */
	private static DateFormat dateFormat =
								new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	/**
	 * Decimal format which is logged in logging simulation messages.
	 */
	@SuppressWarnings("unused")
	private static DecimalFormat myformat = new DecimalFormat("0.0");
	
	/**
	 * Used to hold simulation time at the time of logging.
	 */
	public static long anInstant;
	
	/**
	 * Name of the log file.
	 */
	private static String theLogFileName;
	
	/**
     * \brief Open a logFile at a specified directory, and initialises the time
     * recorder.
     * 
     * @param dirName Directory where the log file should be stored.
     */
	public static void openInitialisationErrorFile(String dirName) 
	{
		try 
		{
			theLogFileName = dirName + File.separator + 
												"Initialisation_Error_Log.txt";
			@SuppressWarnings("unused")
			File test = new File(theLogFileName);
			
			log = new FileOutputStream(theLogFileName);
		} catch (Exception e) 
		{
			System.out.println("Failed to create a log file : "+e);
		}
	}
	
	/**
     * \brief Static method to add message to the log file, even when quiet
     * mode is enabled.
     * 
     * @param message The message that should be appended to the log file.
     */
	public static void writeLog(String message)
	{
		try
		{
			System.out.println(message);
			log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
			log.write((" : "+message+"\n").getBytes());
		}
		catch (Exception e)
		{
			System.out.println("Failed to write into the log file : "+e);
		}
	}
	
	/**
     * \brief Static method to add a specific error message to the log file.
     * 
     * @param exception	Exception associated with this error.
     * @param origin String noting the origin of this error
     */
	public static void writeError(String message,String origin)
	{
		try {
			System.out.println(message);
			log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
			log.write((" Error met in "+origin).getBytes());
			log.write((" : "+message+"\n").getBytes());
		}
		catch (Exception e)
		{
			System.out.println("Failed to write into the log file : "+e);
		}
	}	

	/**
	 * \brief Closes the log file
	 * 
	 * Closes the log file
	 */
	public static void closeFile() 
	{
		try 
		{
			log.close();
		} 
		catch (Exception e) {
		}
	}
	
	/**
	 * \brief Deletes the log file if no errors have been logged
	 * 
	 * Deletes the log file if no errors have been logged
	 */
	public static void deleteFile()
	{
		try
		{
			File f = new File(theLogFileName);
			f.delete();
		}
		catch(Exception e){
		}
	}
}
