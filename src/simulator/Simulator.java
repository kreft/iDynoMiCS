/**
 * \package simulator
 * \brief Package of classes that create a simulator object and capture simulation time.
 * 
 * Package of classes that create a simulator object and capture simulation time. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator;

import idyno.Idynomics;
import idyno.SimTimer;
import java.io.File;
import java.io.ObjectInputStream;
import java.util.*;
import org.jdom.Element;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import de.schlichtherle.io.FileInputStream;
import utils.ExtraMath;
import utils.MTRandom;
import utils.XMLParser;
import utils.ResultFile;
import utils.LogFile;
import utils.Chart;
import povray.PovRayWriter;
import simulator.diffusionSolver.*;
import simulator.geometry.*;
import simulator.reaction.*;
import simulator.agent.*;


/**
 * \brief Top-level class of the simulation core. Used to create and run a simulation
 * 
 * The Simulator class is called by the iDynomics class and creates an object that simulates the conditions specified in a given 
 * protocol file
 * 
 * @since June 2006
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class Simulator 
{

	/**
	 * An XML parser of the protocol file that specifies the conditions under which this simulation is being run
	 */
	public transient XMLParser    _protocolFile;
	
	/**
	 * Boolean that notes that a file defining the agent conditions at the start of the simulation has been specified. Such files are 
	 * normally taken from an evolved state of a similar simulation
	 */
	private boolean               useAgentFile;
	
	/**
	 * Where the parameter useAgentFile is set to true, this must be a parser of an XML file that describes the initial state in which the 
	 * agents should start the simulation
	 */
	private XMLParser             agentFile;
	
	/**
	 * Boolean that notes that a file defining the bulk conditions at the start of the simulation has been specified. Such files are 
	 * normally taken from an evolved state of a similar simulation
	 */
	private boolean               useBulkFile;
	
	/**
	 * Where the parameter useBulkFile is set to true, this must be a parser of an XML file that describes the initial state in which the 
	 * bulk should start the simulation
	 */
	private XMLParser             bulkFile;
	
	/**
	 * Output period at which results are written to file. Specified in the XML protocol file
	 */
	private double                _outputPeriod;
	
	/**
	 * Time counter used to generate a summary report of the previous iteration
	 */
	private double	_lastOutput;
	
	/**
	 * Writer for pov-ray output files. Used in biofilm simulations but not for modelling chemostats
	 */
	public transient PovRayWriter povRayWriter;
	
	/**
	 * Array of result files produced by this simulation
	 */
	public transient ResultFile[] result;
	
	/**
	 * Path to where results files should be stored. Specified in the protocol file
	 */
	private String                _resultPath;

	/**
	 * Boolean that runs the simulation in chemostat conditions
	 */
	public static boolean isChemostat =false;
	
	/**
	 * Boolean stating whether we are using a fluctuating environment. Defaults to false
	 */
	public static boolean isFluctEnv=false;
	
	/**
	 * Flag that triggers the use of the MultiEpiBac and MultiEpisome classes. Added by Sonia Martins Jan 2011
	 */
	public static boolean multiEpi=false;
	
	/**
	 * Flag that notes whether the simulation is 3D or 2D
	 */
	public boolean is3D = false;

	/**
	 * Invasion/Competition simulation. Set to true if the simulation should stop once there is only one species left in the system
	 */
	public static boolean invComp = false;
	
	/**
	 * Boolean controlling simulation stop point - set to false when a LocatedGroup is removed from above the threshold, thus ending 
	 * the simulation. Added by Rob 5/3/11. To date only implemented in DS_SolGrad
	 */
	public boolean continueRunning = true;

	/** 
	 * Timer of the simulation 
	 */
	public static SimTimer        simTimer;

	/** 
	 * Allows user to define a smaller timestep for agent behaviors and interactions. Should always be EQUAL or LOWER than the 
	 * global time step 
	 */
	public double  agentTimeStep;
	
	/**
	 * String that specifies the method of attachment used in the simulation. Can be 'onetime' - where the cells are placed randomly on the 
	 * substratum, or self-attach, where the bacteria moves from the boundary layer to the substratum using a random walk. Defaults to onetime 
	 * as this has been the method employed in previous versions of iDynomics. Added by Kieran, 170513.
	 */
	public String attachmentMechanism = "onetime";

	/**
	 * Specification of all the geometry, bulks, and computational domains in the system being modelled
	 */
	public World                  world;
	
	/**
	 * Dictionary of all the species objects in the specified protocol file. Allows to know the index of an object before its creation
	 */
	public ArrayList<String> 	speciesDic;
	
	/**
	 * Dictionary of all the reaction objects in the specified protocol file. Allows to know the index of an object before its creation
	 */
	public ArrayList<String> 	reactionDic;
	
	/**
	 * Dictionary of all the solver objects in the specified protocol file. Allows to know the index of an object before its creation
	 */
	public ArrayList<String> 	solverDic;
	
	/**
	 * Dictionary of all the particle objects in the specified protocol file. Allows to know the index of an object before its creation
	 */
	public ArrayList<String>	particleDic;
	
	/**
	 * Dictionary of all the solute objects in the specified protocl file. Allows to know the index of an object before its creation
	 */
	public ArrayList<String>	soluteDic;

	
	/**
	 * Array of solute grids - one for each solute specified in the simulation protocol file
	 */
	public SoluteGrid[]           soluteList;
	
	/**
	 * Array of Reaction objects - one for each specified in the simulation protocol file
	 */
	public Reaction[]             reactionList;
	
	/**
	 * Array of Solver objects - one for each solver specified in the simulation protocol file
	 */
	public DiffusionSolver[]      solverList;
	
	/**
	 * Array of Species objects - one for each specifies in the simulation protocol file
	 */
	public ArrayList<Species>     speciesList = new ArrayList<Species>();
	
	/**
	 * List of plasmid carriage information, if creating a Multi Episome simulation. Sonia 07-05-09
	 */
	public ArrayList<String>     plasmidList = new ArrayList<String>();
	
	/**
	 * List of all the scan speeds of all plasmids, if creating a Multi Episome simulation. Sonia 11-10-10
	 */
	public ArrayList<Double>    scanSpeedList = new ArrayList<Double>();

	/** 
	 * Grid where all located agents are stored 
	 */
	public AgentContainer         agentGrid;

	/**
	 * Chart object used by the createCharts method. Deprecated (KA 09/04/13) 
	 */
	private Chart _graphics;
	
	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/
	
	/**
	 * \brief Simulation Constructor. Initialise simulation with protocol file and relevant result and log file paths
	 * 
	 * Creates the simulator object under which a simulation will be run for the conditions stated in the protocol file. This is  
	 * called within the initSimulation method of the Idynomics class. The path to a protocol file and results folder is provided 
	 * and the simulation initialised. Required species, solutes, world, and grids are created and set to their respective initial 
	 * conditions, and result file objects created that will produce statistics at a specified frequency
	 * 
	 * @param protocolFile	Path to the protocol file specifying the conditions under which the simulation is run
	 * @param resultPath	Path to a relevant directory where the simulation results will be stored
	 */
	public Simulator(String protocolFile, String resultPath) 
	{
		try 
		{
			LogFile.chronoMessageIn("System initialisation:");

			// Create pointers to protocolFiles (scenario and agents)
			_protocolFile = new XMLParser(protocolFile);
			_resultPath = resultPath+File.separator;
			
			// Now detect whether initial state files are being used that describe the states the agent and bulk should start the
			// simulation in. If these have been specified, paths to these files will be stored in agentFile and bulkFile respectively
			detectInputFiles(protocolFile);

			// Read in the parameters in the SIMULATOR mark-up of the XML file.
			// This call also deals with setting the time steps (adaptive or specified)
			createSimulator();
			
			// Read in the parameters in the WORLD markup of the XML file. Creates the bulk and computation domains
			createWorld();
			
			// Iterate through the specified solutes in this simulation, creating grids for each and initial concentrations
			createSolutes();
			
			// Now create all the reactions specified in the REACTION mark-up
			createReactions();
			
			// Now to initialise the solvers, specified in the SOLVER mark-up
			createSolvers();
			 
			// Finally, from the protocol file, initialise all the specified species, establishing the required populations of each
			createSpecies();
			
			// KA 09/04/13 - Assume this is a deprecated method as this has always been commented out. Left though for time-being
			//createCharts();

			// bvm 27.1.2009 moved this from above (just below createWorld) to allow better povray outputs
			createFiles(resultPath);
			
			//LogFile.chronoMessageOut("done");

			// Now Describe initial conditions
			if (!Simulator.isChemostat)
			{
				povRayWriter.write(SimTimer.getCurrentIter());
			}	
			
			writeReport();
			 
			// Simulation initialisation complete	 
			 

		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Simulator.CreateSystem(): error met: "+e);
			e.printStackTrace();
			System.exit(-1);
		}

	}

	/**
	 * \brief Method that starts the simulation, calling each timestep until the simulation time is over or the biofilm fills the domain
	 * 
	 * Method that starts the simulation, calling each timestep until the simulation time is over or the biofilm fills the domain
	 */
	public void run() 
	{
		// Rob 5/3/11: added condition belowThreshold. This will be set false whenever
		// a LocatedGroup is removed from above the threshold (if there is one). So
		// far only implemented in DS_SolGrad
		// Rob 7/6/11: changed to a more general continueRunning. Now also (potentially)
		// used to stop a simulation if one species is washed out: see  
		while (simulationIsRunning() && continueRunning)
			step();
		if (!continueRunning) writeReport();
		
		
	}

	/**
	 * \brief Utility method that determines whether the simulation is running
	 * 
	 * Utility method that determines whether the simulation is running
	 * 
	 * @return Boolean stating if the simulation is running (true) or has stopped (false)
	 */
	public boolean simulationIsRunning() 
	{
		return !SimTimer.simIsFinished();
	}

	

	/**
	 * \brief Perform a full iteration of the simulation for the conditions stated in the protocol file
	 * 
	 * Perform a full iteration of the simulation initialised under the conditions in a given protocol file. Accessed via the simulators 
	 * run method
	 */
	public void step() 
	{
		try
		{
			long startTime = System.currentTimeMillis();
	
			// Increment system time
			SimTimer.applyTimeStep();
	
			// Check if new agents should be created
			checkAgentBirth();
			
			
			// Perform diffusion-reaction relaxation
			
			LogFile.chronoMessageIn();
	
			for (DiffusionSolver aSolver : solverList)
				aSolver.initAndSolve();
			LogFile.chronoMessageOut("Solving Diffusion-reaction");
	
	
			//sonia: 25-08-09
			if(isFluctEnv){
				FluctEnv.setEnvCycle(FluctEnv.envNameList.indexOf(FluctEnv.envStatus));
			}
	
			// Perform agent stepping
			agentGrid.step(this);
			LogFile.chronoMessageOut("Simulating agents");
			
			
	
	
			// output result files
			// this will output if we're close to the output period but haven't
			// quite hit it (which happens sometimes due to floating point issues)
			// (this will also output for sure on the last step) 
			if ( ((SimTimer.getCurrentTime()-_lastOutput)
					>= (_outputPeriod - 0.01*SimTimer.getCurrentTimeStep())) ||
					SimTimer.simIsFinished() ) {
				writeReport();
				
				//sonia 26.04.2010
				//only remove the agents from the system after recording all the information about active
				//and death/removed biomass
				agentGrid.removeAllDead();
			}	
			
			
	
			
	
			// Rob 11/2/2012
			// If this is an invComp simulation (default is false), stop if there are fewer than
			// two species remaining
			if (invComp) {
				int nSpecies = speciesDic.size();
				int specAlive = 0;
				
				for (int iSpecies = 0; iSpecies<nSpecies; iSpecies++) {
					if (speciesList.get(iSpecies).getPopulation() > 0)
						specAlive++; 
				}
				if (specAlive < 2)
					continueRunning = false;
			}
	
			SimTimer.updateTimeStep(world);
			LogFile.writeEndOfStep(System.currentTimeMillis()-startTime);
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	
	/**
	 * \brief Processes parameters specified in the SIMULATOR protocol file mark-up, sets the random number generator and creates the data dictionary
	 * 
	 * Within the XML file, the SIMULATOR mark-up defines the parameters that control the overall simulation. These are read in from 
	 * the XML file by this method, initialising the simulation. These parameters include specifying whether this run is under chemostat 
	 * conditions, the setting of the simulation timestep, and time interval between the generation of output files. The random number 
	 * generator is also initialised, either from scratch or from a previous state if this is a re-run. Thirdly, the simulation timer 
	 * is initialised. Finally, a data dictionary is created for use by other utility methods within the simulation
	 */
	public void createSimulator() 
	{
		System.out.print("\t Simulator: ");
		XMLParser localRoot = new XMLParser(_protocolFile.getChildElement("simulator"));

		// Read the flag from protocol file to decide if this is a chemostat run (Added by Sonia)
		isChemostat = localRoot.getParamBool("chemostat");
		
		isFluctEnv = localRoot.getParamBool("isFluctEnv");
		
		multiEpi = localRoot.getParamBool("ismultiEpi");
		
		// invComp - Invasion/Competition simulation - true if the simulation should stop once there is only one species left
		invComp = localRoot.getParamBool("invComp");

		// Read in the agentTimeStep
		agentTimeStep = localRoot.getParamTime("agentTimeStep");
		
		// Read in the method of attachment of agents - whether substratum or boundary layer based.
		// Put in an IF to check as if using protocol files for previous versions, this tag may not be present
		// We want this to remain at the default "onetime" rather than be set to null if the tag is not present
		// KA 170513
		if(localRoot.getParam("attachment")!=null)
		{
			attachmentMechanism = localRoot.getParam("attachment");
		}

		// Now we need to determine if a random.state file exists. This is used to initialise the random number generator
		//sonia 05.2011 bug fix: the code was not finding the random.state file because no path was 
		//being given to the File class to check if the file existed. 
		File randomFile = new File(Idynomics.currentPath+File.separator+"random.state");

		if(randomFile.exists()) 
		{
			/* if a file called random.state exists, the random number generator is initialised using this file.
			 *  this ensures that the random number stream does not overlap when running repeated simulations with the same protocol file
			 *  Chinmay 11/08/2009
			 */				

			FileInputStream randomFileInputStream;
			ObjectInputStream randomObjectInputStream;
			try 
			{
				randomFileInputStream = new FileInputStream(Idynomics.currentPath+File.separator+"random.state");
				randomObjectInputStream = new ObjectInputStream(randomFileInputStream);
				ExtraMath.random = (MTRandom) randomObjectInputStream.readObject();
				LogFile.writeLog("Read in random number generator");
			}
			catch (Exception e) 
			{
				LogFile.writeLog("Simulator.createSimulator() : error met while reading in random number state file" + e);
				System.exit(-1);
			}
			finally
			{
				try
				{
					System.out.println("Random number generator test: "+ExtraMath.random.nextInt());
				}
				catch(java.lang.NullPointerException npe)
				{
					LogFile.writeLog("Random number generator test failed!");
					LogFile.writeLog("See Simulator.createSimulator()");
				}
			}
		}
		else 
		{
			// No random state file exists - initialise the random state generator without any previous information
			long randomSeed;
			try
			{
				String rSeed = localRoot.getParam("randomSeed");
				randomSeed = (long) Integer.parseInt(rSeed);
				LogFile.writeLog("Using random seed in protocol file: "+randomSeed);
			}
			catch (Exception e)
			{
				randomSeed = (long) (Math.random()*Calendar.getInstance().getTimeInMillis());
				LogFile.writeLog("No valid random seed given in protocol file.");
				LogFile.writeLog("Using a randomly generated seed : "+randomSeed);
			}
			
			ExtraMath.random = new MTRandom(randomSeed);
			System.out.println("Random number generator test: "+ExtraMath.random.nextInt());
		}

		// Now to initialise the simulation timer - the timesteps specified in the SIMULATOR markup of the XML file are processed
		// by that timer object
		simTimer = new SimTimer(localRoot);

		// need to reset the time & iterate if we're restarting a run
		if (localRoot.getParamBool("restartPreviousRun")) 
		{

			simTimer.setTimerState(_resultPath+File.separator
					+"lastIter"+File.separator
					+"env_Sum(last).xml");
		}

		// Create dictionary - creates lists of the names of each of the species, particles, solutes, reactions, and solvers in this
		// simulation for reference later
		createDictionary();

		System.out.println("done");
	}

	/**
	 * \brief Create simulation result and POV-Ray files within the required output folder
	 * 
	 * This method creates the required env_State, env_Sum, agent_State, agent_Sum and other required result files within a specified 
	 * result directory. This also initialises a POV-Ray writing object that will produce data that can be represented using POV-Ray
	 *  
	 * @param resultPath	String stating the file path where these files should be created
	 */
	public void createFiles(String resultPath) 
	{
		XMLParser localRoot = new XMLParser(_protocolFile.getChildElement("simulator"));

		_outputPeriod = localRoot.getParamTime("outputPeriod");
		
		// Initialise data files
		// bvm 26.1.2009: added passing of current iterate to output files to 
		// make restarting more robust
		result = new ResultFile[6];
		int currentIter = SimTimer.getCurrentIter();
		result[0] = new ResultFile(resultPath, "env_State", currentIter);
		result[1] = new ResultFile(resultPath, "env_Sum", currentIter);
		result[2] = new ResultFile(resultPath, "agent_State", currentIter);
		result[3] = new ResultFile(resultPath, "agent_Sum", currentIter);

		//sonia 26.04.2010
		//result files for death/removed biomass
		result[4] = new ResultFile(resultPath, "agent_StateDeath", currentIter);
		result[5] = new ResultFile(resultPath, "agent_SumDeath", currentIter);
		
		// Initialise povray files
		// Rob: no need in a chemostat
		if (!Simulator.isChemostat)
		{
			povRayWriter = new PovRayWriter();
			povRayWriter.initPovRay(this, resultPath);
		}
	}

	/**
	 * \brief Builds a set of dictionaries capturing the names of the solutes, particles, reactions, species, and solvers in this specific simulation
	 * 
	 * Creates a set of dictionaries populated with the names of each of the solutes, particles, reactions, species, and solvers, that
	 * were specified in the XML protocol file
	 * 
	 */
	public void createDictionary() 
	{
		LinkedList<Element> list;

		// Build the list of "solutes" markup
		list = _protocolFile.buildSetMarkUp("solute");
		soluteDic = new ArrayList<String>(list.size());
		soluteList = new SoluteGrid[list.size()];
		for (Element aChild : list)
			soluteDic.add(aChild.getAttributeValue("name"));

		// Build the dictionary of "particles"
		list = _protocolFile.buildSetMarkUp("particle");
		particleDic = new ArrayList<String>(list.size());
		for (Element aChild : list)
			particleDic.add(aChild.getAttributeValue("name"));

		
		// Trick to guarantee that the EPS compartment (capsule) is in
		// last position if it exists
		if (particleDic.remove("capsule")) particleDic.add("capsule");

		// Build the dictionary of "reactions"
		list = _protocolFile.buildSetMarkUp("reaction");
		reactionDic = new ArrayList<String>(list.size());
		reactionList = new Reaction[list.size()];
		for (Element aChild : list)
			reactionDic.add(aChild.getAttributeValue("name"));

		// Build the dictionary of "species"
		list = _protocolFile.buildSetMarkUp("species");
		speciesDic = new ArrayList<String>(list.size());
		for (Element aChild : list)
			speciesDic.add(aChild.getAttributeValue("name"));

		// Build the dictionary of "solvers"
		list = _protocolFile.buildSetMarkUp("solver");
		solverDic = new ArrayList<String>(list.size());
		solverList = new DiffusionSolver[list.size()];
		for (Element aChild : list)
			solverDic.add(aChild.getAttributeValue("name"));

	}

	/**
	 * \brief Processes the WORLD mark-up from the protocol file, initialising the bulk and computational domains
	 * 
	 * The world mark-up in the protocol file collects the description of all bulks and computational domains defined in the 
	 * simulation. The bulk mark-up defines a bulk solute compartment that is a source or sink for solutes involved in biofilm growth. 
	 * The computationDomain mark-up defines the spatial region the biofilm will grow in. Only one world may be defined, but this world 
	 * may contain several bulk compartments and computationDomain domains, each with a different name. Though when simulating a 
	 * chemostat scenario, the name of the bulk MUST be chemostat, regardless of the corresponding computationalDomain name. This 
	 * method creates the world properties (e.g. system size and boundary conditions). In terms of boundary conditions, the current 
	 * release allows you to create zero-flow, cyclic, constant and dilution boundaries. To create your own, create a new class 
	 * extending the abstract class "BoundaryCondition"
	 * 
	 */
	public void createWorld() 
	{
		// Creation of the world
		try 
		{
			System.out.print("\t World: \n");
			
			// Get the world markup from the protocol file
			XMLParser parser = new XMLParser(_protocolFile.getChildElement("world"));
			
			// Initialise a simulation world with the setup specified in the protocol file
			world = new World();
			world.init(this, parser);

			// If a initial bulk state is required and set using parameter useBulkFile, we need to recreate these conditions
			if (useBulkFile) recreateBulkConditions();

			System.out.println("\t done");
		} 
		catch (Exception e) 
		{
			// Log the error creating the world
			LogFile.writeLog("Simulator.createWorld() : "+e);
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/**
	 * 
	 * \brief Recreates bulk conditions if an initial state file was specified
	 * 
	 * Get the bulk concentrations from the input file and assign them to the current bulks
	 * 
	 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
	 */
	public void recreateBulkConditions() throws Exception 
	{
		String bulkName;
		int soluteIndex;
		String soluteName;
		XMLParser simulationRoot = new XMLParser(bulkFile.getChildElement("simulation"));

		for (Element aBulkMarkUp : simulationRoot.buildSetMarkUp("bulk")) 
		{

			XMLParser aBulkRoot = new XMLParser(aBulkMarkUp);
			bulkName = aBulkRoot.getAttributeStr("name");

			// check to make sure the bulk exists
			if (!world.containsBulk(bulkName))
				throw new Exception("Bulk "+bulkName+" is not specified in protocol file");

			LogFile.writeLog("\t\tInitializing bulk '"+bulkName+"' from input file.");

			Bulk thisBulk = world.getBulk(bulkName);

			// now set the solutes within this bulk
			for (Element aSoluteMarkUp : aBulkRoot.buildSetMarkUp("solute")) 
			{
				XMLParser aSoluteRoot = new XMLParser(aSoluteMarkUp);

				soluteName = aSoluteRoot.getAttributeStr("name");
				soluteIndex = getSoluteIndex(soluteName);

				// check consistency with protocol file declarations
				if (!soluteDic.contains(soluteName))
					throw new Exception("Solute "+soluteName+" is not in protocol file");

				// finally set the value

				thisBulk.setValue(soluteIndex, Double.parseDouble(aSoluteRoot.getElement().getValue()));
				LogFile.writeLog("\t\tsolute "+soluteName+" is now: "+thisBulk.getValue(soluteIndex));
			}

			LogFile.writeLog("\t\tInitialized bulk '"+bulkName+"' from input file.");
		}

	}

	/**
	 * \brief Initialises each solute specified in the protocol file, creating a solute grid for each and setting the initial concentration
	 * 
	 * In iDynoMiCS each solute exists in its own spatial grid. This method iterates through each specified solute and creates the
	 * necessary grid space, in either 2D or 3D as specified in the protocol file. Once the grid is initialised, the grid is populated 
	 * with values representing the initial concentration levels of the solute, taken from the Sbulk specification in the protocol file
	 * 
	 */
	public void createSolutes() 
	{
		System.out.print("\t Solutes: \n");
		try 
		{
			// Count of solutes and reference to solute list
			int iSolute = 0;
			
			// First get a linked list of all the Solute tags in the XML protocol file
			for (Element aSoluteMarkUp : _protocolFile.buildSetMarkUp("solute")) 
			{
				// Now to parse each solute in turn
				XMLParser aSoluteRoot = new XMLParser(aSoluteMarkUp);
				soluteList[iSolute] = new SoluteGrid(this, aSoluteRoot);
				LogFile.writeLog("\t\t"+soluteList[iSolute].getName()+" ("+soluteList[iSolute].soluteIndex+")");
				iSolute++;
			}
			System.out.println("\t done");
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Simulator.createSolutes() : error met " + e);

			System.exit(-1);
		}
	}

	/**
	 * \brief Create all reactions described in the protocol file
	 * 
	 * This method creates all the reactions that were described in the XML protocol file, storing these in an array of reaction 
	 * objects.
	 */
	public void createReactions() throws Exception 
	{
		Reaction aReaction;

		LogFile.writeLog("\t Reactions: \n");

		// Reference to reactionList array of reactions
		int iReaction = 0;
		
		// Now to go through each reaction in the file
		for (Element aReactionMarkUp : _protocolFile.buildSetMarkUp("reaction")) 
		{
			XMLParser aReactionRoot = new XMLParser(aReactionMarkUp);
			// Create a reaction object
			aReaction = (Reaction) aReactionRoot.instanceCreator("simulator.reaction");
			// Initialise that reation as specified in the XML file
			aReaction.init(this, aReactionRoot);

			// register the created object into the reactions container
			reactionList[iReaction] = aReaction;
			iReaction++;
			
			aReaction.register(aReaction.reactionName, this);
			
			LogFile.writeLog("\t\t"+aReaction.reactionName+" ("+aReaction.reactionIndex+")");
		}
		
		LogFile.writeLog("\t done");
	}

	/**
	 * \brief Creates the pressure solvers specified in the protocol file
	 * 
	 * This method processes the pressure solvers that are specified in the XML protocol file. These are used to compute the steady 
	 * state solute profile within the computational domain
	 */
	public void createSolvers() throws Exception 
	{
		System.out.print("\t Solvers: \n");
		try
		{
			// Now for each specified solver in the XML file
			for (Element aSolverMarkUp : _protocolFile.buildSetMarkUp("solver")) 
			{
				// Initialise the XML parser
				XMLParser parser = new XMLParser(aSolverMarkUp);
	
				// Create the solver,initialise it and register it
				DiffusionSolver aSolver =
					(DiffusionSolver) parser.instanceCreator("simulator.diffusionSolver");
				
				// Create the solver object
				aSolver.init(this, parser);
				
				// Add this to the simulation array of solver objects
				aSolver.register();
				
				System.out.println("\t\t"+aSolver.solverName+" ("+aSolver.solverIndex+")");
			}

			System.out.println("\t done");
		}
		catch(Exception e)
		{
			LogFile.writeLog("Error in Simulator.creatSolvers(): "+e);
		}
	}

	/**
	 * \brief Creates the species, agent grids, and populations required for this simulation.
	 * 
	 * This method is quite complex and is in three distinct parts. The first creates the species and progenitors, and registers these 
	 * in the simulation. The second deals with the creation of the grids in which the species agents exist. Finally, populations of 
	 * each species are created and initialised
	 * 
	 */
	public void createSpecies() throws Exception 
	{
		// THIS PROCESS IS CONDUCTED IN THREE DISTINCT STAGES
		// STAGE 1: CREATE THE SPECIES (AND THE PROGENITOR) AND REGISTER IT
		try 
		{ 
			System.out.print("\t Species: \n");

			XMLParser speciesDefaults;
			// Read in the 'Species Defaults' from the parameter file
			if(_protocolFile.getChildElement("speciesDefaults") != null)
			{
			    speciesDefaults = new XMLParser(_protocolFile.getChildElement("speciesDefaults"));
			}
			else
			{
				speciesDefaults = new XMLParser(new Element("speciesDefaults"));
			}
			
			
			// Now iterate through all species specified in the protocol file
			for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) 
			{
				// Create a new species object for this specification
				Species aSpecies = new Species(this, new XMLParser(aSpeciesMarkUp),speciesDefaults); 
				// Add to the list of species in this simulation
				speciesList.add(aSpecies);
				
				LogFile.writeLog("\t\t"+aSpecies.speciesName+" ("+aSpecies.speciesIndex+")");
			}
		
			System.out.print("\t done\n");
			
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Error in Simulator.createSpecies() first stage");
            e.printStackTrace();
		}
		
		
		// STAGE 2: CREATE THE AGENT GRID
		try 
		{ 
			System.out.print("\t Agent Grid: \n");
			
			// Get the agent grid markup from the XML file
			XMLParser parser = new XMLParser(_protocolFile.getChildElement("agentGrid"));
			
			// Create the grid using specification in the XML file
			agentGrid = new AgentContainer(this, parser, agentTimeStep);
			is3D = agentGrid.is3D;
			System.out.print("\t done\n");

			// Finalise the initialisation of the progenitor
			System.out.print("\t Species progenitor: \n");
			
			for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) 
			{
				parser = new XMLParser(aSpeciesMarkUp);
				
				// Create the progenitor
				getSpecies(parser.getAttribute("name")).getProgenitor().initFromProtocolFile(this,
						parser);
				
				//sonia: creating a list with the plasmid names which will be used afterwards to write the agentSum report
				// Whether MultiEpiBac and MultiEpisome agents are being used 
				if(multiEpi)
				{
					if(parser.getAttribute("class").equals("MultiEpisome"))
					{
						String plName = parser.getAttribute("name");
						Double scanSpeed = parser.getParamDbl("scanSpeed");
						plasmidList.add(plName);
						scanSpeedList.add(scanSpeed);

					}
				}
			}
			
			System.out.print("\t done\n");
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Error in Simulator.createSpecies() second stage");
			e.printStackTrace();
		}
		
		// STAGE 3: CREATION OF POPULATIONS OF THE SPECIES
		// NOTE WE ARE ONLY GOING TO DO STAGE 3 (AT THIS POINT IN DEVELOPMENT) IF ONE TIME ATTACHMENT IS ON
		// FOR SELF ATTACHMENT, ENTRY OF AGENTS BEGINS WHEN THE SIMULATION BEGINS
		
		if(this.attachmentMechanism.equals("onetime"))
		{
			System.out.print("\t Species populations: \n");

			// If there is a previous state file, recreate the species that were present
			if (useAgentFile)
			{
				recreateSpecies();
						
				//sonia: 20-07-09 
				//I've added the line of code below, so that it is possible to create new agents from the protocol file when restarting
				//from a previous simulation (environment) using a new protocol file
				checkAgentBirth();
	
			}
			else
			{
				checkAgentBirth();
			}
		
			System.out.print("\t done\n");
		}
			
	}

	/**
	 * \brief Create agents from a file describing individually all the agents species by species. Used if useAgentFile is true
	 * 
	 * It is possible to initialise the simulation with information in an agent file. If this is the case, the species will need to be 
	 * created from information in this file. This method performs a simulation initialisation in that case
	 * 
	 * @throws Exception	Exception thrown if this file cannot be read correctly or file is not consistent with protocol file
	 */
	public void recreateSpecies() throws Exception 
	{
		int spIndex, counterSpecies;
		SpecialisedAgent progenitor;
		
		// Initialise a parser of the XML agent file
		XMLParser simulationRoot = new XMLParser(agentFile.getChildElement("simulation"));
		counterSpecies = 0;
		
		for (Element aSpeciesMarkUp : simulationRoot.buildSetMarkUp("species")) 
		{
			XMLParser aSpeciesRoot = new XMLParser(aSpeciesMarkUp);
			spIndex = getSpeciesIndex(aSpeciesRoot.getAttribute("name"));
			
			// check consistency with protocol file declarations
			boolean isConsistent = (speciesList.get(counterSpecies).speciesName.equals(aSpeciesRoot
					.getAttributeStr("name")));
			
			if (!isConsistent) throw new Exception("Agent input file is inconsistent with protocol file: ");
			
			// Else process agents description
			String dataSource = aSpeciesRoot.getElement().getContent(0).toString();
			String[] allAgentData = null;
			try 
			{
				allAgentData = dataSource.split(";\n");
				// this removes the '[Text: \n' line from the first item in the list
				// so that all agents will be treated equally
			} 
			catch (Exception e)
			{
				LogFile.writeLog("Simulator.recreateSpecies() : problem splitting up data");
			}
			
			allAgentData[0] = allAgentData[0].substring(8);
			
			progenitor = speciesList.get(spIndex).getProgenitor();
			
			// don't use the last index because it is a string with only ']'
			for (int i = 0; i < allAgentData.length-1; i++)
				(progenitor.sendNewAgent()).initFromResultFile(this, allAgentData[i].split(","));
			
			LogFile.writeLog(speciesList.get(counterSpecies).speciesName+" : "
					+speciesList.get(counterSpecies).getPopulation()
					+" agents created from input file.");
			counterSpecies++;
		}
	}

	/**
	 * \brief Introduces a new agent to the simulation, initialises the specified birthday and initial area parameters
	 * 
	 * Introduces a new agent to the simulation, initialises the specified birthday and initial area parameters
	 */
	public void checkAgentBirth() 
	{	
		XMLParser parser;
		int spIndex;
		boolean creatingAgents = false;

		// Now go through each species
		for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) 
		{
			parser = new XMLParser(aSpeciesMarkUp);
			spIndex = getSpeciesIndex(parser.getAttribute("name"));
			
			// KA May 2013: We are now adding additional methods of attachment. Thus the old code has been moved from here into its
			// own method. A switch is used to determine which method of attachment is being employed.
			// Although there are only 2 at the moment (in version 1.2), a switch is useful if further methods are added later
			
			if(this.attachmentMechanism.equals("onetime"))
			{
				creatingAgents = this.oneTimeAttachmentAgentBirth(parser,spIndex);
			}
			else if(this.attachmentMechanism.equals("selfattach"))
			{
				creatingAgents = this.selfAttachmentAgentBirth(parser,spIndex);
			}
			
		}
		if (creatingAgents) agentGrid.relaxGrid();
	}

	/**
	 * \brief Creates a specified number of agents within a specified area of the computation domain, where one time attachment is being employed
	 * 
	 * From Version 1.2 of iDynoMiCS, different methods of attachment are being introduced alongside the original one-time adhesion to the 
	 * substratum. This method captures the case in the original version of the tool, where a number of cells are created at a random 
	 * location within a set region of the domain.
	 * 
	 * @param parser	The XML tags contained within the "species" markup
	 * @param spIndex	The index of this species in the simulation dictionary
	 * @return	A boolean noting whether any agents have been created
	 */
	public boolean oneTimeAttachmentAgentBirth(XMLParser parser, int spIndex)
	{
		boolean creatingAgents = false;
		
		for (Element aInitMarkUp : parser.buildSetMarkUp("initArea")) 
		{
			parser = new XMLParser(aInitMarkUp);
			
			if (SimTimer.isDuringNextStep(parser.getParamTime("birthday"))) 
			{
				this.agentGrid.agentIter = agentGrid.agentList.listIterator();
				speciesList.get(spIndex).createPop(parser);
				creatingAgents = true;
			}
		}
		
		return creatingAgents;
	}
	
	/**
	 * \brief For self-attachment scenarios, creates the required number of agents per timestep. These migrate to the substratum or biofilm surface
	 * 
	 * For self attachment scenarios, cells are initialised on the boundary layer surface rather than the substratum. iDynoMiCS then models 
	 * the agents move towards the biofilm surface using a run and tumble motion. If the cell moves towards the biofilm or surface, the 
	 * move continues until the cell attaches. If the move takes the cell back towards the bulk, we assume the cell is back in the mix and 
	 * introduce a further cell for a new attempt at attaching.
	 * 
	 * @param parser	The XML tags contained within the "species" markup
	 * @param spIndex	The index of this species in the simulation dictionary
	 * @return	A boolean noting whether any agents have been created
	 */
	public boolean selfAttachmentAgentBirth(XMLParser parser, int spIndex)
	{
		boolean creatingAgents = false;
		
		// Now check the injection period. 
		if(simTimer.getCurrentTime()>=speciesList.get(spIndex).cellInjectionStartHour
			&& simTimer.getCurrentTime()<=speciesList.get(spIndex).cellInjectionEndHour)
			{
				// We're in a cell injection period. In this time, we model cells being injected into the domain. These may then
				// attach at a different frequency to when this injection is off
			creatingAgents = 
					create_Required_Number_Of_SelfAttaching_Agents(speciesList.get(spIndex).injectionOnAttachmentFrequency, 
																	parser,spIndex);
			}
		else
		{
			// Injection is off, attach at the different rate
			creatingAgents = 
					create_Required_Number_Of_SelfAttaching_Agents(speciesList.get(spIndex).injectionOffAttachmentFrequency, 
																	parser,spIndex);
		}	
		
		return creatingAgents;
	}
	
	/**
	 * \brief Creates the required number of cells at the boundary layer in self-attachment scenarios. Differs dependent on whether within an injection period
	 * 
	 * For self-attachment scenarios, it is possible to model injection of cells during the simulation, as well as periods that follow 
	 * where no cells are injected. In these cases, the number of cells injected from the boundary layer will differ (to model the fact 
	 * that the attachment rate may be different in these periods). This rate is supplied to this function, and the required number of 
	 * cells created
	 * 
	 * @param cellAttachmentFrequency	The number of cells to create at the boundary layer per hour, for this injection period
	 * @param parser	The XML object for this species. From the simulation protocol file
	 * @param spIndex	The index of this species in the simulation species list
	 * @return	A boolean noting whether any agents were created in this step
	 */
	public boolean create_Required_Number_Of_SelfAttaching_Agents(double cellAttachmentFrequency, XMLParser parser, int spIndex)
	{
		boolean creatingAgents = false;
		
		// To keep things consistent, the global timestep is specified in hours, and the number of cells injected is also specified 
		// in hours. Thus we need to work out we are injecting the correct number of cells here.
		
		int wholeAgentsThisTimeStep = (int)(cellAttachmentFrequency * SimTimer.getCurrentTimeStep());
		double remainder = (cellAttachmentFrequency * SimTimer.getCurrentTimeStep()) - (double)wholeAgentsThisTimeStep;
		
		if(remainder>0)
		{
			// increment the remainder for this species
			speciesList.get(spIndex).newAgentCounter = speciesList.get(spIndex).newAgentCounter+remainder;
		}
		
		// Now see if the remainder has totalled a new cell
		if(speciesList.get(spIndex).newAgentCounter >= 1)
		{
			wholeAgentsThisTimeStep++;
			speciesList.get(spIndex).newAgentCounter--;
		}
		
		// Now create these agents
		if(wholeAgentsThisTimeStep>0)
		{
			this.agentGrid.agentIter = agentGrid.agentList.listIterator();
			speciesList.get(spIndex).createBoundaryLayerPop(parser,wholeAgentsThisTimeStep);
			creatingAgents = true;
			
		}
		
		return creatingAgents;
	}

	/**
	 * \brief Creates reports of concentrations of solutes and biomass on the default gird and a description of all agents
	 * 
	 * Generate a report with concentrations of solutes and biomass on the default grid and a report with the exhaustive description 
	 * of all agents. Each report is a new file with a new index
	 */
	public void writeReport() 
	{
		// Update saving counters and file index
		_lastOutput = SimTimer.getCurrentTime();
		int currentIter = SimTimer.getCurrentIter(); // bvm added 26.1.2009

		// first restart log file to avoid non-write trouble
		LogFile.reopenFile();

		try 
		{
			/* Grids and environment ______________________________________ */
			// env_State
			result[0].openFile(currentIter);
			// env_Sum
			result[1].openFile(currentIter);
			// bvm added 16.12.08

			//sonia:chemostat
			if(Simulator.isChemostat)
			{
				//sonia:chemostat
				//I've modified refreshBiofilmGrids()
				soluteList[0].getDomain().refreshBioFilmGrids();
			}
			else
			{
				// output the biofilm thickness data
				//sonia 12.10.09				

				Double [] intvals;
				StringBuffer value = new StringBuffer();
				for (Domain aDomain : world.domainList) 
				{
					aDomain.refreshBioFilmGrids();
					intvals = aDomain.getInterface();

					value.append("<thickness domain=\""+aDomain.domainName+"\" unit=\"um\">\n");
					value.append("\t<mean>"+(ExtraMath.mean(intvals))+"</mean>\n");
					value.append("\t<stddev>"+(ExtraMath.stddev(intvals, false))+"</stddev>\n");
					value.append("\t<max>"+(ExtraMath.max(intvals))+"</max>\n");
					value.append("</thickness>\n");

				}

				result[0].write(value.toString());
				result[1].write(value.toString());

			}

			// Add description of each solute grid
			for (SoluteGrid aSG : soluteList) 
			{
				aSG.writeReport(result[0], result[1]);
			}

			// Add description of each reaction grid
			for (Reaction aReac : reactionList) 
			{
				aReac.writeReport(result[0], result[1]);
				
				// KA - August 2013 - For each reaction, calculate the production/uptake of each solute
				// Thus for each output period, we can report an estimate of production/uptake of each solute
				// Where outputPeriod > simulation step, the output remains the last simulation step. The user can interpolate
				// from this if required.
				aReac.calculateSoluteChange();
			}
			

			// Add description of each species grid
			agentGrid.writeGrids(this, result[0], result[1]);

			// Add description of total biomass
			for (Domain aDomain : world.domainList) 
			{
				aDomain.refreshBioFilmGrids();
				aDomain.getBiomass().writeReport(result[0], result[1]);
			}

			// KA AUGUST 2013
			// Now we're going to add to ENV_STATE the amount of solute produced or consumed in this domain in this step
			summariseSoluteProductionOrUptake();			
			
			// Add description of bulks
			// THIS COMPLETES THE ENV_STATE SOLUTE COUNTS, UNDER THE BULK NAME TAGS
			for (Bulk aBulk : world.bulkList) 
			{
				aBulk.writeReport(result[1]);
			}

			// Close EnvState and EnvSum
			result[0].closeFile();
			result[1].closeFile();

		} 
		catch (Exception e) 
		{
			LogFile.writeError("System description of grids failed:"+e.getMessage(),
			"Simulator.writeReport()");
		}

		try 
		{
			/* Agents ____________________________________________________ */
			result[2].openFile(currentIter);
			result[3].openFile(currentIter);
			result[4].openFile(currentIter);
			result[5].openFile(currentIter);

			agentGrid.writeReport(this, result[2], result[3]);
			agentGrid.writeReportDeath(this, result[4], result[5]);

			result[2].closeFile();
			result[3].closeFile();
			result[4].closeFile();
			result[5].closeFile();

			// Rob 15/2/2011: No need to write povray if it's a chemostat
			if (!Simulator.isChemostat) povRayWriter.write(currentIter);

			LogFile.writeLog("System description finalized");

		} 
		catch (Exception e) 
		{
			LogFile.writeError("System description of agents failed:"+e.getMessage(),
			"Simulator.writeReport()");
		}

	}
	
	/**
	 * \brief Calculates the production or uptake of each solute over all reactions for this timestep and writes them to the env_state file
	 * 
	 * Calculates the production or uptake of each solute over all reactions for this timestep and writes them to the env_state file. This 
	 * gives some indication of the solute change in this simulation. Note that unlike other aspects of this simulation, the output here 
	 * remains the solute figures from just the previous step should the output period be greater than the simulation timestep. The user 
	 * will have to interpolate from there if this is the case
	 */
	public void summariseSoluteProductionOrUptake()
	{
		// KA AUGUST 2013
		// Now we're going to add to ENV_STATE the amount of solute produced or consumed in this domain, in this step
		StringBuffer value = new StringBuffer();
		value.append("<globalProductionRate>\n");
		// Need to iterate through each solute to do this
		for(int iSolute=0;iSolute<soluteDic.size();iSolute++)
		{
			value.append("<solute name=\"").append(soluteDic.get(iSolute));
			value.append("\" unit=\"").append("g.hour");
			value.append("\">\n");
			
			// Now we need the total production and uptake for this solute across all the reactions.
			// Each reaction holds this info
			double soluteReactionFigure = 0;
			for (Reaction aReac : reactionList) 
			{
				soluteReactionFigure+=aReac._soluteProductionOrUptake[iSolute];
			}

			// Now we have the value, write it to the buffer
			value.append(soluteReactionFigure);
			
			value.append("</solute>");
		}
		value.append("</globalProductionRate>\n");
		// Append the gloabl solute values to the env state file
		result[1].write(value.toString());
	}
	
	/**
	 * \brief Deprecated internal chart creation method. Assumed deprecated as never called in code (KA 09/04/13)
	 * 
	 * Internal chart creation method. Assumed deprecated as never called in code (KA 09/04/13)
	 */
	public void createCharts() 
	{

		_graphics = new Chart("Simulation outputs");
		_graphics.setPath(_resultPath);

		XYSeriesCollection[] graphSet = new XYSeriesCollection[2];
		String[] xLegend = {"Time(h)","Time(h)"};
		String[] yLegend =  {"Conc (g.L-1)","Conc (g.L-1)"};

		// Create a chart for Solutes		
		graphSet[0]=new XYSeriesCollection();
		graphSet[1]=new XYSeriesCollection();
		int nSolute = soluteDic.size();
		int nSpecies = speciesDic.size();
		for (int iSolute = 0; iSolute<nSolute; iSolute++)
			graphSet[0].addSeries(new XYSeries(soluteDic.get(iSolute)));		
		for (int iSpecies = 0; iSpecies<nSpecies; iSpecies++)		
			graphSet[1].addSeries(new XYSeries(speciesDic.get(iSpecies)));	

		_graphics.init(graphSet,xLegend,yLegend);

		_graphics.pack();
		_graphics.setVisible(true);

	}
	
	/**
	 * \brief Deprecated method. Not used in this version of iDynoMiCS (KA 09/04/13)
	 * 
	 * Deprecated method. Not used in this version of iDynoMiCS (KA 09/04/13)
	 */
	public void updateChart() 
	{
		int nSolute = soluteDic.size();
		int nSpecies = speciesDic.size();

		for (int iSolute = 0; iSolute<nSolute; iSolute++)
			_graphics.updateChart(0,iSolute, SimTimer.getCurrentTime(), world.getBulk("tank")
					.getValue(iSolute));

		for (int iSpecies = 0; iSpecies<nSpecies; iSpecies++)
			_graphics.updateChart(1,iSpecies, SimTimer.getCurrentTime(),
					speciesList.get(iSpecies).getPopulation());

		_graphics.repaintAndSave();


	}


	/**
	 * \brief Find a Species on the basis of its specified string name
	 * 
	 * Returns a species on the basis of its specified string name
	 * 
	 * @param aSpeciesName
	 * @return the speciesIndex
	 */
	public int getSpeciesIndex(String aSpeciesName) 
	{
		return speciesDic.indexOf(aSpeciesName);
	}

	/**
	 * \brief Returns a species object for the given species name
	 * 
	 * Returns a species object for the given species name
	 * 
	 * @param aSpeciesName	Text string noting the name of the species object required
	 * @return	The species object of that specified name
	 */
	public Species getSpecies(String aSpeciesName) 
	{
		return speciesList.get(getSpeciesIndex(aSpeciesName));
	}

	/**
	 * \brief Return the integer reference specifying the position of a given solute in the solute dictionary
	 * 
	 * The createDictionary method created a dictionary of all solutes used in the simulation. This method is used to return the 
	 * index of a particular solute, found by the given solute name
	 * 
	 * @param aSoluteName	The name of the solute whos position in the dictionary is being sought
	 * @return	Integer specifying the position in the dictionary / list
	 */
	public int getSoluteIndex(String aSoluteName) 
	{
		return soluteDic.indexOf(aSoluteName);
	}

	public SoluteGrid getSolute(String aSoluteName) {
		return soluteList[getSoluteIndex(aSoluteName)];
	}

	public int getReactionIndex(String aReactionName) {
		return reactionDic.indexOf(aReactionName);
	}

	/**
	 * \brief Returns a reaction object for a reaction of a given name
	 * 
	 * Returns a reaction object for a reaction of a given name
	 * 
	 * @param aReactionName	The name of the reaction to return
	 * @return	A reaction object that describes a specified reaction
	 */
	public Reaction getReaction(String aReactionName) 
	{
		return reactionList[getReactionIndex(aReactionName)];
	}

	/**
	 * \brief Returns an integer reference to the position of a solver in the simulation solver dictionary
	 * 
	 * @param aSolverName	The name of the solver to reference
	 * @return	An integer value specifying the position of this solver in the simulation dictionary
	 */
	public int getSolverIndex(String aSolverName) 
	{
		return solverDic.indexOf(aSolverName);
	}

	/**
	 * \brief Returns a solver object of a specified name
	 * 
	 * Returns a solver object of a specified name
	 * 
	 * @param aSolverName	Name of the solver to return from the array of solvers
	 * @return	Solver object of the specified name
	 */
	public DiffusionSolver getSolver(String aSolverName) 
	{
		int solInd = getSolverIndex(aSolverName);
		if (solInd >= 0)
			return solverList[solInd];
		else
			return null;
	}

	public int getParticleIndex(String particleName) {
		return particleDic.indexOf(particleName);
	}

	/**
	 * \brief	Determines whether a previous run is being restarted, reading in the respective last state if this is the case 
	 * 
	 * Check whether we will initialize from agent and bulk input files. If the "restartPreviousRun" param is true, this will also set
	 * the correct files for reading in the last state
	 * 
	 * @param protocolFile	The name of the protocol file for which the simulation is being restarted
	 */
	public void detectInputFiles(String protocolFile) {

		// first check whether we are restarting from a previous run
		XMLParser restartInfo = new XMLParser(_protocolFile.getChildElement("simulator"));

		if (restartInfo.getParamBool("restartPreviousRun")) {
			// if this is true, then we set the input files as the last files
			// that were output

			useAgentFile = true;
			useBulkFile = true;

			agentFile = new XMLParser(_resultPath+File.separator
					+"lastIter"+File.separator
					+"agent_State(last).xml");
			bulkFile = new XMLParser(_resultPath+File.separator
					+"lastIter"+File.separator
					+"env_Sum(last).xml");

			LogFile.writeLog("Restarting run from previous state in directory: "
					+_resultPath);

			return;
		} 

		// otherwise just do things as usual, but only if input is specified
		if (_protocolFile.getChildElement("input")==null) return;

		XMLParser input = new XMLParser(_protocolFile.getChildElement("input"));

		useAgentFile = input.getParamBool("useAgentFile");
		if (useAgentFile) {
			String agentFileName = input.getParam("inputAgentFileURL");
			// construct the input file name using the path of the protocol file
			int index = protocolFile.lastIndexOf(File.separator);
			agentFileName = protocolFile.subSequence(0, index+1)+agentFileName;

			agentFile = new XMLParser(agentFileName);
			LogFile.writeLog("Using agent input file: "+agentFileName);
		}

		useBulkFile = input.getParamBool("useBulkFile");
		if (useBulkFile) {
			String bulkFileName = input.getParam("inputBulkFileURL");
			// construct the input file name using the path of the protocol file
			int index = protocolFile.lastIndexOf(File.separator);
			bulkFileName = protocolFile.subSequence(0, index+1)+bulkFileName;

			bulkFile = new XMLParser(bulkFileName);
			LogFile.writeLog("Using bulk input file: "+bulkFileName);
		}

	}

	/**
	 * \brief Return the path to the simulation result files
	 * 
	 * Return the path to the simulation result files
	 * 
	 * @return	String value stating the path to the result files
	 */
	public String getResultPath() {
		return _resultPath;
	}

}