/**
 * \package idyno
 * \brief Package of classes used to launch iDynomics
 * 
 * Package of classes used to launch and iDynoMiCS simulation, and to update the package to the latest stable release. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the 
 * following URL  "http://www.cecill.info".
 */
package idyno;

import java.util.LinkedList;

import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import de.schlichtherle.io.FileOutputStream;
import simulator.Simulator;
import soft_update.Check_Release_Status;
import utils.InitialisationErrorLog;
import utils.XMLParser;
import utils.ExtraMath;
import utils.LogFile;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.ObjectOutputStream;
import utils.ResultFile;

/** \brief Main class to run to launch the iDynoMiCS tool
*
* Launch iDynomics using this class. Syntax: To launch a single simulation idyno.Idynomics
* \\protocol\\protocolFileURL - To launch a batch of simulations idyno.Idynomics
* \\protocol\\FilesDirectory
* 
* Please make sure you are aware of the licensing agreement when you are running the iDynomics tool. This software is governed by 
* the CeCILL license under French law and  abiding by the rules of distribution of free software.  You can  use, modify and/ or 
* redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  
* "http://www.cecill.info". 
* 
* @since June 2006
* @version 1.2
* @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de)
* @author Laurent Lardon (lardonl@supagro.inra.fr)
* @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
* @author João Xavier (xavierj@mskcc.org)
* @author Cristian Picioreanu (C.Picioreanu@tudelft.nl)
* @author Jan-Ulrich Kreft (j.kreft@bham.ac.uk)
* @author Rob Clegg (rjc096@bham.ac.uk)
* @author Kieran Alden (k.j.alden@bham.ac.uk)
* 
**/
public class Idynomics 
{
	
	/**
	 * Version number of this iteration of iDynoMiCS - required by update procedure
	 */
	public static double version_number  = 1.2;
	
	/**
	 * Frame to display information regarding a software update, if there is one available
	 */
	private static JFrame updateFrame;
	
	/**
	 * Text panel within the update frame that takes text from a URL on the iDynomics server (history.html) to display what 
	 * changes have been made in this update
	 */
	private static JEditorPane infoPane;
	
	/**
	 * Scroll pane to contain the description of the changes should these be numerous
	 */
	private static JScrollPane scp;
	/**
	 * Ok button on the panel to start the update
	 */
	private static JButton ok;
	/**
	 * Cancel button on the updater panel to run without installing the update
	 */
	private static JButton cancel;
	
	/**
	 * Panel for laying out the updater screen
	 */
	private static JPanel pan1;
	
	/**
	 * Panel for laying out the updater screen
	 */
	private static JPanel pan2;
	
	/**
	 * The program arguments sent in with this run (which have yet to be processed) - passed round such that these can be processed
	 * after the update check
	 */
	public static String[] programArgs;
	
	/**
	 * Store new version number (if there is one) globally as required by action handlers
	 */
	public static double newVersionNumber;
	
	/**
	 * An array of the protocol XML file names to be processed, each describing a particular simulation scenario
	 */
	private static String[] _protocolFile;

	/**
	 * An array of the paths to each protocol file to be processed
	 */
	private static String[] _activePath;

	/**
	 * A double used to calculate how long the simulation has taken to run under a set of specified conditions
	 */
	public static double    begin;
	
	/**
	 * A simulation object that runs the simulator under the conditions specified in the protocol file
	 */
	public static Simulator aSimulator;
	
	
	/**
	 * When this is true, only the absolute essentials are displayed and written to the log file
	 */
	public static boolean quietMode = false;

	/**
	 * parameter set such that the simulator knows the current path to read/write the random number file
	 */
	public static String currentPath;
	
	/**
	 * Log to hold any issues with simulation initialisation. Created as initialisation happens prior to creation of any log file
	 */
	public static InitialisationErrorLog initErrorLog;

	/**
	 * \brief - Main method used to start iDynoMiCS
	 * 
	 * This method starts iDynoMiCS and launches the required methods to process the protocol files provided by the user. Firstly, 
	 * these command line or Eclipse dialog box arguments are processed, such that the tool is aware of the protocol files it needs 
	 * to process. With this complete, the tool performs a check that tests whether it is the latest stable release of iDynoMiCS. If this 
	 * is the case, the simulation starts as normal. If not, the user has the opportunity to update iDynoMiCS and relaunch the tool. Note 
	 * that this check can be turned off in the protocol file, and thus a further check is included in this file to determine if the
	 * check is enabled. 
	 * 
	 * @param args Simulation arguments passed from the command line
	 */
	public static void main(String[] args) 
	{
		
		try
		{
			initErrorLog = new InitialisationErrorLog();
			
			processArguments(args);
			
			// Now set the current path to that of the first protocol file to be processed (activePath is set by processArguments)
			currentPath = _activePath[0];
			System.out.println(currentPath);
			// Open the initialisation error log
			
			InitialisationErrorLog.openInitialisationErrorFile(currentPath);
			
			// The first thing we're going to do is check for an iDynoMiCS stable release update - if the user has specified this.
			// The user can determine whether this is done or not in the protocol files they are submitting to iDynoMiCS.
			
			// KA - for sake of saving checking time and avoiding problems with bulk runs, the decision has been made to only check 
			// the FIRST protocol file. Although this could be deemed a bit risky, as the user may select a check in the other protocol 
			// files, this is important, as checking every file will (a) take time, and (b) could cause the run of simulations to be 
			// halted if the first few all say don't check for updates, yet one protocol file then says to. So we should be instructing 
			// users that if they do wish to check for updates, this needs to be noted in the first protocol file that they specify. 
			// The below method checks the protocol file, and if update check is enabled, performs the check as to whether this is the 
			// latest version. If not, this returns true. If it is the latest version or the update check is set to false, 
			// false is returned
			
			boolean updateCheck = checkForReleaseUpdate(_activePath[0],_protocolFile[0]);
			
			if(updateCheck)
			{
				// Call the method to show the update box and offer the user the chance to update
				show_Update_Box(args, Check_Release_Status.getWhatsNew());	
			}
			else
			{
				// Start processing the protocol files
				// Note this has been included as a call rather than the code be placed here, as the start must be accessible from the
				// update screen too - as the user may choose not to update.
				startIDynomics();
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}

	}
	

	/** \brief Method to show that an iDynomics update is available
	 * 
	 * This method creates a dialog box on the screen that shows that an update is available for the iDynomics package. 
	 * This dialog box will contain information from a file history.html, on the iDynomics web server, that details the changes 
	 * included within this new version. The user is given two options - to either install the update (after which they will need to 
	 * restart iDynomics), or proceed with the simulation without installing the update
	 * 
	 * @param args	The program arguments sent in to iDynomics
	 * @param updateInfoMessage	The list of changes in this update, from history.html
	 * 
	 */
	public static void show_Update_Box(String[] args, String updateInfoMessage)
	{
		// Store the program arguments so these can be processed if an update is not performed
		programArgs = args;
		
		// Create the update frame
		updateFrame = new JFrame();
		updateFrame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
		updateFrame.setTitle("New Update Found");
		pan1 = new JPanel();
		pan1.setLayout(new BorderLayout());
		 
		pan2 = new JPanel();
		pan2.setLayout(new FlowLayout());
		 
		infoPane = new JEditorPane();
		infoPane.setContentType("text/html");
		infoPane.setText(updateInfoMessage);
		 
		scp = new JScrollPane();
		scp.setViewportView(infoPane);
		 
		// Create the Update button, and set this to start the jar that performs the update. Note that the jar must be in the
		// updater folder
		
		ok = new JButton("Update");
		ok.addActionListener( new ActionListener(){
		 
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				// Run the updater
				String[] run = {"java","-jar","src/lib/Update_IDynomics.jar"};
				try 
				{
					Runtime.getRuntime().exec(run);
			    } 
				catch (Exception ex) 
				{
					ex.printStackTrace();
				}
				// Exit once this is started, as the updater may not be able to update some of the running classes if this is still
				// running
				System.exit(0);		
			}
		});
		 
		// Add a cancel button incase the user does not want to run the update - if this is the case iDynomics will start
		cancel = new JButton("Run Without Update");
		cancel.addActionListener( new ActionListener(){
		 
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				updateFrame.dispose();
				try
				{
					Idynomics.startIDynomics();
				}
				catch(Exception ex)
				{
					ex.printStackTrace();
				}
		    }
		});
		
		// Finish setting up the panels
		pan2.add(ok);
		pan2.add(cancel);
		pan1.add(pan2, BorderLayout.SOUTH);
		pan1.add(scp, BorderLayout.CENTER);
		updateFrame.add(pan1);
		updateFrame.pack();
		updateFrame.setVisible(true);
		updateFrame.setSize(500, 500);
	}
	
	
	/**
	 * \brief Launches a simulation of each of the protocol files provided
	 * 
	 * This method iterates through each of the protocol files provided and launches a simulation under each of these conditions.  
	 * Once this is complete, the random number state file for that simulation is also written to file. Note that from version 1.2 this 
	 * method is new - but all this code previously existed in the main method. It has been moved here to work with the updater - the user 
	 * may get the update screen but choose not to update. If this is the case, they need to be located somewhere where they can launch 
	 * the simulation. This was only possible by creating a new method 
	 */
	public static void startIDynomics() throws Exception
	{
		System.out.println("Starting iDynomics");
		
		// Start the simulation
		// Now to process each protocol file in turn
		for (int iSimul = 0; iSimul<_protocolFile.length; iSimul++)
		{
			if (initSimulation(iSimul)) 
			{
				// If open, close the initialisation error log
				try
				{
					InitialisationErrorLog.closeFile();
					InitialisationErrorLog.deleteFile();
					
				}
				catch(Exception e)
				{
					// No need to perform an action
				}
				
				launchSimulation(iSimul);
				
				
						
				/* The following lines write out the random number state file at the end of each simulation
				 * The reason this is done here, is because it is guaranteed to be the absolute last thing that's done before
				 * the next simulation is called. Chinmay 11/8/2009 */
				
				try 
				{
					FileOutputStream randomFileOutputStream = new FileOutputStream(currentPath+File.separator+"random.state");
					ObjectOutputStream randomObjectOutputStream = new ObjectOutputStream(randomFileOutputStream);
					randomObjectOutputStream.writeObject(ExtraMath.random);
					randomObjectOutputStream.close();
					randomFileOutputStream.close();
					LogFile.writeLogAlways("Wrote random number generator");
					LogFile.closeFile();
				}
				catch(Exception e) 
				{
					LogFile.writeLogAlways("Idynomics.main(): error met while writing out random number state file" + e);
					LogFile.closeFile();
					System.exit(-1);					
				}
				
				
				
				
			}
		}
	}

	/**
	 * \brief Method to select protocol files from a file selection dialog if run in Eclipse
	 * 
	 * When run from the command line (using the RunIDyno.py file), the simulation expects XML protocol file names to be provided.
	 * If however the simulation is run from Eclipse (which is possible, especially for those developing it), the user would have to
	 * either specify the file name as a Run Argument, or provide this during RunTime. This method allows the latter, presenting the 
	 * simulation user with a selection box from which to choose one or more protocol files on which simulations should be run. These 
	 * are returned within a linked list for later processing
	 * 
	 * @return a linked list, protocolFiles, of XML files selected from the dialog box
	 */
	public static LinkedList<File> chooseFile() 
	{
		/**
		 * protocolFiles - a linked list that will contain the files selected from the dialog box
		 */
		LinkedList<File> protocolFiles = new LinkedList<File>();

		// Open a FileChooser window in the current directory
		JFileChooser chooser = new JFileChooser(""+System.getProperty("user.dir")+"/protocol");
		chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		
		// Allow the user to select multiple files
		chooser.setMultiSelectionEnabled(true);

		if (chooser.showOpenDialog(null)==JFileChooser.APPROVE_OPTION) 
		{
			// Go through each file and add to the protocolFiles list
			// This is performed using the utility method listDirectory
			for (java.io.File aFile : chooser.getSelectedFiles()) 
			{
				listDirectory(aFile, "xml", protocolFiles);
			}
		}

		return protocolFiles;
	}

	/**
	 * \brief Processes the command line or Eclipse protocol file selection arguments, forming a list of simulations to be run
	 * 
	 * iDynoMiCS can be run in a couple of different ways. Yet both require the user to specify one or a set of XML protocol files that 
	 * describe the simulation that will be run. This method is used in both Eclipse and command line version of iDynoMiCS, and processes 
	 * this input. For a command line run (from RunIDyno.py for example), the user will specify the protocol files on the command line. For 
	 * an eclipse run, these will be specified via a file selection dialog box. Either way, these are then processed into two string arrays, 
	 * one containing the paths to the protocol files to be run and another containing the names of the protocol files. Once this is complete. 
	 * iDynoMiCS can be started such that each is processed in turn.
	 * 
	 * @param args	The protocol file names specified on the command line (not used in an eclipse run)
	 */
	static public void processArguments(String[] args) 
	{
		LinkedList<File> allFiles;

		// We're going to build a list of the protocol files upon which simulations will be run.
		// As there are two ways of running the simulation (from the RunIDyno script where a protocol file must be provided,
		// and eclipse where a protocol file can be selected from a window) we need to check whether a protocol file (or files) have
		// been provided as arguments
		
		switch (args.length) 
		{
			case 0:
				// The simulation is being run in Eclipse
				// A protocol file (or files) need to be specified - and this is done using a file selection window
				allFiles = chooseFile();
				break;
			default:
				// The assumption here is that protocol files have been provided as command line arguments. RunIDyno.py checks that
				// this is the case before the simulation is run
				
				// Initialise the linked list (if wondering why not done in above case, this was processed in chooseFile)
				allFiles = new LinkedList<File>();
				
				// Now go through the arguments, and process each using the listDirectory utility method (could be a directory or file)
				for (int iFile = 0; iFile<args.length; iFile++) 
				{
					listDirectory(new File(args[iFile]), "xml", allFiles);
				}
		}

		// Now to remove any agent_State and env_Sum files from the list
		// KA 26/03/13 - When I looked at this I saw no requirement to do this (as the files were never there, as the XML filter 
		// made sure of this). I have however left this here for the moment - this should be reviewed in the future
		
		// Get the number of files in the list
		int nfile = allFiles.size();
		// Keep track of the files we've checked
		int nchecked = 0;
		// While there are still files to check
		while (nchecked < nfile) 
		{
			// Remove the file if contains agent_State or env_Sum
			if (allFiles.get(nchecked).getName().contains("agent_State") ||
					allFiles.get(nchecked).getName().contains("env_Sum")) 
			{
				allFiles.remove(nchecked);
				nfile = allFiles.size();
			} 
			else 
			{
				// Just inc to show the file at this reference in the list has been checked
				nchecked = nchecked + 1;
			}
		}

		// Now to build the list of path and protocol files that will be processed in this run
		int nProtocol = allFiles.size();
		
		// Initialise the arrays that will store (a) the path to the protocol file and (b) the names of the protocol file
		_activePath = new String[nProtocol];
		_protocolFile = new String[nProtocol];

		// Iterate through each protocol file
		for (int iFile = 0; iFile<nProtocol; iFile++) 
		{
			// Get the full path to this file and add to the array
			_activePath[iFile] = allFiles.get(iFile).getParent()+java.io.File.separator;
			
			//Edd: Added check for null _activePath, which was occuring when the launch scripts
			//were invoked directly from the folder containing the protocol files in use.
			if(_activePath[iFile].equals("null"+java.io.File.separator))
			{
				_activePath[iFile]="."+java.io.File.separator;
			}
			
			// Now add the name of the protocol file to the array of those to be processed
			_protocolFile[iFile] = allFiles.get(iFile).getName();
			
			// Tell the user what we're doing
			System.out.println("Initializing with protocol file: "+_activePath[iFile]
			                                                                   +_protocolFile[iFile]);
		}
	}

	/**
	 * \brief Create the result directories and simulator object for the simulation to be run
	 * 
	 * Called for each protocol file specified in the command line or Eclipse arguments. This method first of all checks whether the 
	 * specified simulation is new or a restart of a previous simulation. If the latter is the case the log file is reopened and a log 
	 * added stating this is the restart of a previous simulation. If a new simulation, the relevant result file directories are created, 
	 * a log file created for population as the simulation progresses, and the protocol file copied into the results folder for ease of 
	 * reference alongside results. Once the required folders are created (where necessary), a simulator object is created for this run. 
	 * Should this be successful, true is returned to note the simulation was initialised ok. If not, the main program is informed such 
	 * that the run is not performed.
	 * 
	 * @param iSimul	The number that references a protocol file in the array of protocol files being processed
	 */
	public static boolean initSimulation(int iSimul) 
	{
		// bvm 23.01.09: added code
		// the restart code works by reading in the "restartPreviousRun" parameter
		// from the protocol file; if true, then the run is continued from the
		// previous state.
		//
		// this approach means a user must edit the output protocol file in the
		// results folder and then choose it (FROM THE RESULT DIRECTORY) when 
		// starting a new run
		// 

		String resultDir;
		String resultFullURL;

		if (checkForRestart(_activePath[iSimul],_protocolFile[iSimul])) 
		{
			resultDir = _activePath[iSimul];
			resultFullURL = resultDir+File.separator+_protocolFile[iSimul];

			// Create a log file
			LogFile.openFile(resultDir);

			LogFile.writeLogAlways("Restarting run from previous results.");
		} 
		else 
		{
			// DEFAULT - NON-RESTARTING METHOD

			// bvm 10.2.2009
			// create the output file directory using the protocol file title too
			// The results are stored in the results directory - thus take the protocol file address and replace protocol with results
			resultDir = _activePath[iSimul].replace("protocol", "results")
			+ _protocolFile[iSimul].replace(".xml", "(")
			+ LogFile.getDateFileName() + ")";

			// Create this directory
			new File(resultDir+File.separator+"lastIter").mkdirs();

			// Copy the protocol file into the result directory
			resultFullURL = resultDir+File.separator+_protocolFile[iSimul];
			ResultFile.copyFile(_activePath[iSimul]+_protocolFile[iSimul], resultFullURL);

			// Create a log file
			LogFile.openFile(resultDir);
		}
		

		// Check whether QuietMode has been specified - if so note this in the log file
		quietMode = checkForQuietMode(_activePath[iSimul]+_protocolFile[iSimul]);
		if (quietMode) 
		{
			LogFile.writeLogAlways("Quiet mode: on");
		}
		
		// Create the simulator
		
		try 
		{
			aSimulator = new Simulator(_activePath[iSimul]+_protocolFile[iSimul], resultDir);
			LogFile.writeLogAlways("Initialization ("+resultFullURL+"):OK");
			return true;
		} 
		catch (Exception e) 
		{
			LogFile.writeLogAlways("Initialization ("+resultFullURL+"):ERROR");
			return false;
		}
		
	}

	/**
	 * \brief Checks the protocol file to see whether a previous simulation is being restarted
	 * 
	 * This method checks the protocol file provided to determine whether a previous run is being restarted. If this is the case, the 
	 * original results file folder and log files are resused
	 * 
	 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
	 * 
	 * @param activePath	Path to the first protocol file being processed
	 * @param protocolFile	Name of the protocol file that is being checked
	 * @return boolean stating whether this is a restarted simulation or not
	 */
	public static boolean checkForRestart(String activePath,String protocolFile) 
	{
		// KA June 2013 - error was occuring here if XML file was incorrect - the log hasn't yet been opened
		// As it is difficult to open the log without the result address, the call has been changed to one that tries to open the XML
		// file but outputs an error message to the screen if this is incorrect
		XMLParser inputFile = new XMLParser(activePath,protocolFile);
		XMLParser restartInfo = new XMLParser(inputFile.getChildElement("simulator"));

		return restartInfo.getParamBool("restartPreviousRun");
	}

	/**
	 * \brief Checks the protocol file to determine whether quiet mode has been enabled
	 * 
	 * This method checks the protocol file provided to determine whether quiet mode has been specified. If this is the case, the 
	 * simulation will note this in the log file
	 * 
	 * @param protocolFile	The path to the protocol file that is being checked
	 * @return	boolean stating whether quiet mode is enabled or not
	 */
	public static boolean checkForQuietMode(String protocolFile) 
	{
		XMLParser inputFile = new XMLParser(protocolFile);
		XMLParser info = new XMLParser(inputFile.getChildElement("simulator"));
		
		return info.getParamBool("quietMode");
	}
	
	/**
	 * \brief Checks the protocol file to determine whether the user wishes to check for iDynoMiCS update, and if so performs the check
	 * 
	 * This method checks the protocol file provided to determine whether the user has turned iDynoMiCS update off. If this is the 
	 * case, the simulation will not check that it is the most up to date version. However, if this is true, the simulation performs a check
	 * against the online repository as to whether this is the latest version. If it is not the latest version, true is returned, and the
	 * user will be shown a screen detailing the changes in the new version, and given the option as to whether to update. If either this is 
	 * the latest version or the user has turned updating off, false will be returned. Note that if checkForReleaseUpdate is not in the
	 * protocol file, FALSE is returned and no check will be performed
	 * 
	 * @author Kieran Alden (k.j.alden@bham.ac.uk)
	 * 
	 * @param activePath	Path to the first protocol file being processed
	 * @param protocolFile	The path to the protocol file that is being checked
	 * @return	boolean stating whether an update needs to be installed (true) or whether the software is up to date or update checking is off (false)
	 */
	public static boolean checkForReleaseUpdate(String activePath,String protocolFile) throws Exception
	{
		XMLParser inputFile = new XMLParser(activePath,protocolFile);
		XMLParser info = new XMLParser(inputFile.getChildElement("simulator"));
		
		if(info.getParamBool("checkForReleaseUpdate"))
		{
			// Check whether this is is the latest version
			return iDyno_New_Released_Version_Check();
		}
		else
		{
			return false;
		}
	}
	
	/**
	 * \brief Checks whether this version of iDynoMiCS is the latest stable release
	 * 
	 * This method checks with the online repository whether this version of iDynoMiCS is up to date. If this is not the case, true 
	 * is returned and the user will be prompted to update iDynoMiCS if they wish. If this is the latest version, false will be returned 
	 * and the simulation starts as normal
	 * 
	 * @author Kieran Alden (k.j.alden@bham.ac.uk)
	 * @return	boolean stating whether an update is required (true) or not (false)
	 * @throws Exception	Exception generated if a problem checking release status online (but this is dealt with elsewhere should this occur)
	 */
	public static boolean iDyno_New_Released_Version_Check() throws Exception
	{
		// Determine if there is a new release available using Check_Release_Status
		// As check release status checks for an internet connection, this exception is caught by returning
		// a version number of 999. If this is returned, there is no internet to check for updates, and thus the simulation
		// should start as normal
			
		newVersionNumber = Double.parseDouble(Check_Release_Status.getLatestVersion());
			
		// Now determine if an update needs to be done
		return (newVersionNumber > version_number && newVersionNumber != 999.0);
	}
	
	/**
	 * \brief Executes the simulation object for a particular protocol file
	 * 
	 * Launches a simulation for each of the provided protocol files, and logs how long the simulation has taken to run under those
	 * conditions
	 * 
	 * @param iSimul	An integer reference to the protocol file name array
	 */
	public static void launchSimulation(int iSimul) 
	{
		try 
		{
			// Get the current system time in milliseconds - used to calculate length of run
			begin = System.currentTimeMillis();
			// Start the simulation
			aSimulator.run();
			
			
			begin = Math.round(System.currentTimeMillis()-begin);

			// Calculate how long the simulation took
			String time = ExtraMath.toString(begin/1e3/60, false);
			
			// Log this in the log file
			LogFile.writeLogAlways("Simulation succesfully achieved in "+time+" minutes.");
			

		} 
		catch (Exception e) 
		{
			// Log the error
			System.out.println("At Idynomics:launch simulation error met :" + e);
			LogFile.writeLogAlways("Simulation failed. " + e);
		}
	}

	/**
	 * \brief Processes the files and directories selected from the dialog box when iDynoMiCS is run in Eclipse
	 * 
	 * When iDynoMiCS is run, the user has to specify protocol files. This is done either via a dialog box in Eclipse, or from the 
	 * command line if run from RunIDyno.py. This method takes the files and directories that have been specified either way and processes these, 
	 * adding the relevant files to a linked list. There are a few considerations here (a) If only one file is selected, the method 
	 * checks this is the correct file type (specified by input argument filter) and if so adds this to the list. (b) If the selection 
	 * is a directory, the whole directory is added to the list, again after a pre-check that the files are of the correct file format
	 * 
	 * @param aDirectory	The file(s) or directory selected from the dialog box or specified on the command line for processing
	 * @param filter	The file format of the protocol files (basically the filter to apply to these files)
	 * @param fileList	The linked list that the protocol files should be added to
	 */
	public static void listDirectory(java.io.File aDirectory, String filter, LinkedList<File> fileList) 
	{
		// First check the selection if not a directory
		// This IF is called if multiple files are selected as well
		if (!aDirectory.isDirectory() && aDirectory.getName().contains("."+filter)) 
		{
			// if it was just a file, add the file as the only member in the list
			fileList.addLast(aDirectory);
		} 
		else 
		{
			// Now we're processing a selected directory - we want to add all files that are XML to the list
			// So first filter the list
			java.io.File[] list = aDirectory.listFiles(new utils.Jfilter(filter));
			
			// Now add each to the linked list
			for (int i = 0; i<list.length; i++) 
			{
				fileList.addLast(list[i]);
			}
		}
	}
	
	/**
	 * \mainpage iDynoMiCS Version 1.2 Code Reference Manual
	 * 
	 * The iDynoMiCS software simulates the growth of microbial communities. iDynoMiCS is written in Java, and uses XML files for input 
	 * and output. Input files allow the users to specify conditions, microbial species, and other parameters without the need to 
	 * become a programmer. Also the positions and other variables of all individual microbes can be specified in another input file. 
	 * This allows one to easily specify many different types of simulations. iDynoMiCS writes plain-text XML files as output, and these 
	 * may be processed using any number of software tools (though we provide some general post-processing routines that run in Matlab 
	 * and are working on R scripts). In addition to XML files, iDynoMiCS also writes scene description input files for POV-Ray, a 
	 * virtual photography software which is used to render 3-D images of the simulated communities.
	 *
	 * iDynoMiCS is under constant development by several different teams working in the realm of microbial ecology. If you are 
	 * interested in contributing to the further development of this software, whether this be through the addition of new 
	 * functionality, new routines for post-processing simulations, or additional PDE solvers or other numerical algorithms, 
	 * please read the information in our <a href="http://www.github.com/kreft">Github repository wiki</a>
	 *
	 * iDynoMiCS has been a collaborative effort between an increasing number of people: Laurent Lardon, Brian Merkey, Joao Xavier, 
	 * Andreas Dötsch, Barth F Smets, Cristian Picioreanu, Rob Clegg, Sonia Martins, Katrin Bohl, Kieran Alden, Jan-Ulrich Kreft 
	 * and others.
	 *   
	 *   <br>
	 *  \section getiDyno Getting iDynoMiCS
	 *   
	 * The tool can be downloaded in two different ways. Should you wish to use the tool to perform experimentation, then download 
	 * the Latest Stable release <a href="http://www.biosciences-labs.bham.ac.uk/kreftlab"> here </a>. If however you wish to contribute to the 
	 * future development of the tool (which we encourage), you can download the most up to date version of the source code from our 
	 * <a href="http://www.github.com/kreft">Github repository</a>. Full documentation that details how to contribute to iDynoMiCS development using 
	 * this repository is available on the iDynoMiCS wiki, also located in the Github repository. 
	 *   
	 *   <br>
	 * \section iDynoMail Join the iDynoMiCS Mailing List
	 * 
	 * While we make the code available under a GNU GPL type open source licence called CeCILL, we kindly request that you join the iDynoMiCS mailing 
	 * lists when you download the code. This will help us find out how many scientists are using iDynoMiCS, which we will use to 
	 * support future grant applications to further develop the code, which in turn is in your interest. We will also use the email 
	 * lists to send updates concerning the ongoing development of iDynoMiCS. <a href="http://www.birmingham.ac.uk/generic/idynomics/mailing-list.aspx">Click here for information on joining these lists</a> 
	 *   
	 * <br>  
	 * \section idyno Help with Using iDynoMiCS
	 * 
	 * The best place to start with iDynoMiCS is to complete the tutorial. This can be found in our  <a href="http://www.github.com/kreft">Github repository</a> 
	 * and a PDF can be found in the 'Tutorial' folder in the iDynoMiCS download.
	 * 
	 * Should there be any issues with using iDynoMiCS, you can seek support in two ways. 
	 * 
	 * The first of these is to ask your question via the <a href="https://github.com/kreft/iDynoMiCS/issues">Issue tracking system</a> 
	 * on our <a href="http://www.github.com/kreft">GitHub repository</a>.  There is more information on how to use the Issue tracking 
	 * system in our wiki part of the GitHub repository. 
	 * 
	 * Secondly, you may send questions to the developers via the mailing list idynomics@lists.bham.ac.uk. The developers will aim to 
	 * respond via email as swiftly as time allows. Your query will also be logged as an Issue in our GitHub repository with the 
	 * appropriate reply, in order that this remains available for other users who may have the same question. 
	 *   
	 *    
	 */
	
	/**
	 * \page page1 Useful Links
	 * 
	 * \section links Useful Links
	 * 
	 * <a href="http://www.idynomics.org">1. The iDynoMiCS Website</a>
	 * 
	 * <a href="https://github.com/kreft">2.. iDynoMiCS GitHub Repository</a>
	 * 
	 * <a href="https://github.com/kreft/iDynoMiCS/wiki">3. iDynoMiCS Wiki</a>
	 * 
	 * <a href="https://github.com/kreft/iDynoMiCS/wiki/iDynomics-Tutorial">4. iDynoMiCS Tutorial</a>
	 * 
	 * <a href="https://github.com/kreft/iDynoMiCS/issues">5. iDynoMiCS Issue Tracking System</a>
	 * 
	 * <a href="http://www.birmingham.ac.uk/generic/idynomics/mailing-list.aspx">6. iDynoMiCS Mailing List</a>
	 * 
	 * <a href="http://www.birmingham.ac.uk/generic/idynomics/projects/index.aspx">7. Projects and Publications</a>
	 * 
	 * <a href="http://www.biosciences-labs.bham.ac.uk/kreftlab">7. Kreft Lab Website</a>
	 * 
	 * 
	 */
	
}