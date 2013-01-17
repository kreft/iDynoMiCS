/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
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



public class LogFile {

	public static FileOutputStream log;
	private static DateFormat      dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static DecimalFormat   myformat   = new DecimalFormat("0.0");
	public static long             anInstant;
	private static String theLogFileName;

	/**
     * Open a logFile and initialise time recorder
     * @param dirName
     */
	public static void openFile(String dirName) {
		try {
			//log = new FileOutputStream(dirName+File.separator+"log.txt");
			
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
		} catch (Exception e) {
			System.out.println("Failed to create a log file");
		}
	}

	/**
     * Static method to add message to the log file, even when quiet mode is enabled
     * @param message
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
     * Static method to add message to the log file, so long as Idynomics.quietMode is false (default)
     * @param message
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
     * Static method to add message to the log file
     * @param message
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

	public static void closeFile() {
		try {
			log.close();
		} catch (Exception e) {
		}
	}
	
	// bvm added 03.09.09
	// this closes and reopens the file in append mode to try and keep
	// logfile updates even when disk location is not writeable
	public static void reopenFile() {
		try {
			closeFile();
			log = new FileOutputStream(theLogFileName, true);
			//System.out.println("Reopened log file.");
		} catch (Exception e) {
		}
	}

	public static void chronoMessageIn(String message) {
		anInstant = System.currentTimeMillis();
		System.out.println(message);
	}

	public static void chronoMessageIn() {
		anInstant = System.currentTimeMillis();
	}

	public static long chronoMessageOut(String message) {
		long value = anInstant;
		anInstant = System.currentTimeMillis();
		value = anInstant-value;

		//System.out.println("\t "+message+" done in "+value/1000+" sec");
		writeLog(message+" done in "+value/1000+" sec");
		return value;
	}

	public static long chronoMessageOut() {
		long value = anInstant;
		anInstant = System.currentTimeMillis();
		return anInstant-value;
	}

	public static void writeEndOfStep(double length) {
		// Simulation started for xxx minutes
		double simTime = (System.currentTimeMillis()-Idynomics.begin)/1000/60;
		if (!Idynomics.quietMode) System.out.println("");

		LogFile.writeLog("Computation time :"+myformat.format(simTime)+" minute(s) \n \t -> Iter "
		        +SimTimer.getCurrentIter()+", Time: "+SimTimer.getCurrentTime()+" achieved in "
		        +myformat.format(length/1000)+" sec \n");
	}
	
	public static String getDateFileName(){
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmm");
		return dateFormat.format(Calendar.getInstance().getTime());		
	}
}
