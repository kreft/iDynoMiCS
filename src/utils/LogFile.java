/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import idyno.Idynomics;
import idyno.SimTimer;

import java.io.File;
import java.io.FileOutputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
// log file where events are recorder
import java.util.Calendar;
import java.text.DecimalFormat;

import Jama.Matrix;


/**
 * \brief Creates and updates a log file for the running simulation
 *
 * Creates and updates a log file for the running simulation
 */
public class LogFile 
{

	/**
	 * Output stream which log messages are written to
	 */
	public static FileOutputStream log;
	
	/**
	 * Format of the date which is used in logging simulation messages
	 */
	private static DateFormat      dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	/**
	 * Decimal format which is logged in logging simulation messages
	 */
	private static DecimalFormat   myformat   = new DecimalFormat("0.0");
	
	/**
	 * Used to hold simulation time at the time of logging
	 */
	public static long             anInstant;
	
	/**
	 * Name of the log file
	 */
	private static String theLogFileName;

	/**
     * \brief Open a logFile at a specified directory, and initialises the time recorder
     * 
     * Open a logFile at a specified directory, and initialises the time recorder
     * 
     * @param dirName	Directory where the log file should be stored
     */
	public static void openFile(String dirName) 
	{
		try 
		{
			// bvm 26.1.2009: create new log files for restarting runs
			int iter = 0;
			File test = new File(dirName+File.separator+"log"+iter+".txt");
			// if the file exists, iterate on the number until the name is free
			while (test.exists()) {
				iter++;
				test = new File(dirName+File.separator+"log"+iter+".txt");
			}
			theLogFileName = dirName+File.separator+"log"+iter+".txt";
			
			log = new FileOutputStream(theLogFileName);
		} catch (Exception e) 
		{
			System.out.println("Failed to create a log file");
		}
	}


	/**
     * \brief Static method to add message to the log file, even when quiet mode is enabled
     * 
     * Static method to add message to the log file, even when quiet mode is enabled
     * 
     * @param message	The message that should be appended to the log file
     */
	public static void writeLogAlways(String message) {
		try {
			System.out.println(message);
			log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
			log.write((" : "+message+"\n").getBytes());
		} catch (Exception e) {
			System.out.println("Failed to write into the log file");
		}
	}
	
	/**
     * \brief Static method to add message to the log file, so long as Idynomics.quietMode is false (default)
     * 
     * Static method to add message to the log file, so long as Idynomics.quietMode is false (default)
     * 
     * @param message	The message that should be appended to the log file
     */
	public static void writeLog(String message) {
		try {
			if (!Idynomics.quietMode){
				System.out.println(message);
				log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
				log.write((" : "+message+"\n").getBytes());
			}
		} catch (Exception e) {
			System.out.println("Failed to write into the log file");
		}
	}
	
	/**
	 * \brief Write a data matrix of a given name into the log file
	 * 
	 * Write a data matrix of a given name into the log file
	 * 
	 * @param name	Name of the matrix being written in the log file
	 * @param M	Data matrix being written into the log file
	 */
	public static void writeMatrix(String name, Matrix M) {
		try {
			LogFile.writeLog("Matrix "+name);
			StringBuilder row = new StringBuilder();
			for (int i = 0; i<M.getRowDimension(); i++){
				for (int j = 0; j<M.getColumnDimension(); j++){
					row.append("     "+M.get(i,j));
				}
				System.out.println(row);
				log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
				log.write((" : "+row+"\n").getBytes());
				row.delete(0, row.length());
			}
		} catch (Exception e) {
			System.out.println("Failed to write into the log file");
		}
	}
	
	/**
     * \brief Static method to add a specific error message to the log file
     * 
     * Static method to add a generated error message to the log file
     * 
     * @param message	Error message that should be appended to the log file
     * @param origin	String noting the origin of this error
     */
	public static void writeError(String message,String origin) {
		try {
			System.out.println(message);
			log.write(dateFormat.format(Calendar.getInstance().getTime()).getBytes());
			log.write((" Error met in "+origin).getBytes());
			log.write((" : "+message+"\n").getBytes());
		} catch (Exception e) {
			System.out.println("Failed to write into the log file");
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
	 * \brief Closes and reopens the log file in append mode to ensure log file updates are kept when disk is not writeable
	 * 
	 * Added by bvm added 03.09.09. Method to close and reopen the log file in append mode, to try and keep logfile updates even 
	 * when disk location is not writeable
	 */
	public static void reopenFile() 
	{
		try 
		{
			closeFile();
			log = new FileOutputStream(theLogFileName, true);
			//System.out.println("Reopened log file.");
		} 
		catch (Exception e) 
		{
		}
	}

	/**
	 * \brief Output a log message with the current time
	 * 
	 * Output a log message with the current time
	 * 
	 * @param message	Message to be appended into the log file
	 */
	public static void chronoMessageIn(String message) {
		anInstant = System.currentTimeMillis();
		System.out.println(message);
	}

	/**
	 * \brief Determine the current time.
	 * 
	 * Determine the current time.
	 */
	public static void chronoMessageIn() {
		anInstant = System.currentTimeMillis();
	}

	/**
	 * \brief Writes a message into the log file with the time taken to perform an action
	 * 
	 * Writes a message into the log file with the time taken to perform an action
	 * 
	 * @param message	The message to be appended into the log file
	 * @return	The amount of time taken to perform an action
	 */
	public static long chronoMessageOut(String message) 
	{
		long value = anInstant;
		anInstant = System.currentTimeMillis();
		value = anInstant-value;

		//System.out.println("\t "+message+" done in "+value/1000+" sec");
		writeLog(message+" done in "+value/1000+" sec");
		return value;
	}

	/**
	 * \brief Determine the amount of time that has passed between two timepoints and return
	 * 
	 * Determine the amount of time that has passed between two timepoints and return
	 * 
	 * @return	Long value noting the milliseconds difference between two timepoints
	 */
	public static long chronoMessageOut() {
		long value = anInstant;
		anInstant = System.currentTimeMillis();
		return anInstant-value;
	}

	/**
	 * \brief Writes log message at the end of the step summarising the computational time taken in this step
	 * 
	 * Writes log message at the end of the step summarising the computational time taken in this step
	 * 
	 * @param length	Length of time in hours that it has taken to perform this step
	 */
	public static void writeEndOfStep(double length) {
		// Simulation started for xxx minutes
		double simTime = (System.currentTimeMillis()-Idynomics.begin)/1000/60;
		if (!Idynomics.quietMode) System.out.println("");

		LogFile.writeLog("Computation time :"+myformat.format(simTime)+" minute(s) \n \t -> Iter "
		        +SimTimer.getCurrentIter()+", Time: "+SimTimer.getCurrentTime()+" achieved in "
		        +myformat.format(length/1000)+" sec \n");
	}
	
	/**
	 * \brief Returns the date and time as a string for use in constructing the log file title
	 * 
	 * Returns the date and time as a string for use in constructing the log file title
	 * 
	 * @return	String containing the year,month,day,underscore,hour, and minute, for appending to the log file name
	 */
	public static String getDateFileName(){
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmm");
		return dateFormat.format(Calendar.getInstance().getTime());		
	}
}
