/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * top-level class of the simulation core. It is used to create and run a
 * simulation; this class is called by the class Idynomics.
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
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

public class Simulator {

	/* ______________ PROPERTIES ____________________________ */

	/** Protocol and optional input XML files describing the simulation * */
	public transient XMLParser    _protocolFile;
	private XMLParser             agentFile;
	private boolean               useAgentFile;
	private XMLParser             bulkFile;
	private boolean               useBulkFile;

	/** Variables used to generate result files and back-ups * */
	private double                _outputPeriod, _lastOutput;
	public transient PovRayWriter povRayWriter;
	public transient ResultFile[] result;
	private String                _resultPath;

	/** Variable used to test the chemostat condition */

	public static boolean isChemostat =false;
	public static boolean isFluctEnv=false;
	//sonia 01/2011
	// flag that triggers the use of the MultiEpiBac and MultiEpisome classes
	public static boolean multiEpi=false;
	public boolean is3D = false;

	public static boolean invComp = false;
	
	// Rob 5/3/11: this will be set to false when a LocatedGroup is removed from above the threshold
	// (if there is one!). So far only implemented in DS_SolGrad.
	public boolean continueRunning = true;

	/** Timer of the simulation */
	public static SimTimer        simTimer;

	/** agentTimeStep */
	// to make it available to AgentContainer class and at the same time be readable from the simulator
	// markup in the protocol file
	public double  agentTimeStep;

	/* The world lists the geometry and the bulks connected to the system */
	public World                  world;

	/*
	 * Dictionary of objects : allows to know the index of an object before its
	 * creation
	 */
	public ArrayList<String>      particleDic, soluteDic, speciesDic, reactionDic, solverDic;

	/* MAIN CONTAINERS ______________________________________________________ */

	/* List of solutes, reactions, solvers and species */
	public SoluteGrid[]           soluteList;
	public Reaction[]             reactionList;
	public DiffusionSolver[]      solverList;
	public ArrayList<Species>     speciesList = new ArrayList<Species>();

	//sonia:07-05-09
	//used to write output information about plasmid carriage
	public ArrayList<String>     plasmidList = new ArrayList<String>();

	//sonia 11.10.2010 list with all the scan speeds of all plasmids
	public ArrayList<Double>    scanSpeedList = new ArrayList<Double>();

	/* Grid where all located agents are stored */
	public AgentContainer         agentGrid;

	private Chart _graphics;

	/* __________________________ CONSTRUCTOR _______________________________ */

	/** Open the protocol_file and initialise the system with */
	public Simulator(String protocolFile, String resultPath) {

		try {
			LogFile.chronoMessageIn("System initialisation:");

			// Create pointers to protocolFiles (scenario and agents)
			_protocolFile = new XMLParser(protocolFile);
			_resultPath = resultPath+File.separator;
			detectInputFiles(protocolFile);

			createSimulator();
			createWorld();
			// createFiles was here
			createSolutes();
			createReactions();
			createSolvers();
			createSpecies();
			//createCharts();

			// bvm 27.1.2009 moved this from above to allow better povray outputs
			createFiles(resultPath);

			LogFile.chronoMessageOut("done");

			// Describe initial conditions
			if (!Simulator.isChemostat)
				povRayWriter.write(SimTimer.getCurrentIter());
			
			writeReport();
			

		} catch (Exception e) {
			LogFile.writeLog("Simulator.CreateSystem(): error met: "+e);
			System.exit(-1);
		}

	}

	/* _____________________ TOP-LEVEL METHODS ______________________________ */
	/**
	 * This is the method that starts the simulation. It is called by the class
	 * with the "main" method
	 */
	public void run() {
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
	 * Check, if the simulation should continue or stop.
	 * 
	 * @return true, if simulation is still running
	 */
	public boolean simulationIsRunning() {
		return !SimTimer.simIsFinished();
	}

	public void createCharts() {

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
	 * Perform a full iteration
	 */
	public void step() {
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
		}	

		//sonia 26.04.2010
		//only remove the agents from the system after recording all the information about active
		//and death/removed biomass
		agentGrid.removeAllDead();

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

	public void updateChart() {
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

	/* ________________ INITIALIZATION PROCEDURE ___________________________ */

	/**
	 * 
	 */
	public void createSimulator() {
		System.out.print("\t Simulator: ");
		XMLParser localRoot = new XMLParser(_protocolFile.getChildElement("simulator"));

		//sonia: read the flag from protocol file to decide if this is a chemostat run
		isChemostat = localRoot.getParamBool("chemostat");
		isFluctEnv = localRoot.getParamBool("isFluctEnv");
		multiEpi = localRoot.getParamBool("ismultiEpi");
		invComp = localRoot.getParamBool("invComp");

		agentTimeStep = localRoot.getParamTime("agentTimeStep");

		//sonia 05.2011 bug fix: the code was not finding the random.state file because no path was 
		//being given to the File class to check if the file existed. 
		File randomFile = new File(Idynomics.currentPath+File.separator+"random.state");

		if(randomFile.exists()) {
			/* if a file called random.state exists, the random number generator is initialised using this file.
			 *  this ensures that the random number stream does not overlap when running repeated simulations with the same protocol file
			 *  Chinmay 11/08/2009
			 */				

			FileInputStream randomFileInputStream;
			ObjectInputStream randomObjectInputStream;
			try {
				randomFileInputStream = new FileInputStream(Idynomics.currentPath+File.separator+"random.state");
				randomObjectInputStream = new ObjectInputStream(randomFileInputStream);
				ExtraMath.random = (MTRandom) randomObjectInputStream.readObject();
				LogFile.writeLog("Read in random number generator");
			}

			catch (Exception e) {
				LogFile.writeLog("Simulator.createSimulator() : error met while reading in random number state file" + e);
				System.exit(-1);
			}
			finally
			{
				try
				{
					System.out.println("rng test: "+ExtraMath.random.nextInt());
				}
				catch(java.lang.NullPointerException npe)
				{
					System.out.println("Random number generator test failed!");
					System.out.println("See Simulator.createSimulator()");
				}
			}
		}
		else {
			//Chinmay - added MTRandom.java to utils and changed the RNG to use that class instead. 11/08/2009
			System.out.println("No random file here!");

			ExtraMath.random = new MTRandom((long) localRoot.getParamDbl("randomSeed"));
		}

		simTimer = new SimTimer(localRoot);

		// need to reset the time & iterate if we're restarting a run
		if (localRoot.getParamBool("restartPreviousRun")) {

			simTimer.setTimerState(_resultPath+File.separator
					+"lastIter"+File.separator
					+"env_Sum(last).xml");
		}

		createDictionary();

		System.out.println("done");
	}

	public void createFiles(String resultPath) {
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
		if (!Simulator.isChemostat){
			povRayWriter = new PovRayWriter();
			povRayWriter.initPovRay(this, resultPath);
		}
	}

	public void createDictionary() {
		LinkedList<Element> list;

		/* Build the list of "solutes" markup _________________________ */
		list = _protocolFile.buildSetMarkUp("solute");
		soluteDic = new ArrayList<String>(list.size());
		soluteList = new SoluteGrid[list.size()];
		for (Element aChild : list)
			soluteDic.add(aChild.getAttributeValue("name"));

		/* Build the dictionary of particles _________________________ */
		list = _protocolFile.buildSetMarkUp("particle");
		particleDic = new ArrayList<String>(list.size());
		for (Element aChild : list)
			particleDic.add(aChild.getAttributeValue("name"));

		// Trick to guarantee that the EPS compartment (capsule) is in
		// last position if it exists
		if (particleDic.remove("capsule")) particleDic.add("capsule");

		/* Build the dictionary of reactions _________________________ */
		list = _protocolFile.buildSetMarkUp("reaction");
		reactionDic = new ArrayList<String>(list.size());
		reactionList = new Reaction[list.size()];
		for (Element aChild : list)
			reactionDic.add(aChild.getAttributeValue("name"));

		/* Build the dictionary of species ___________________________ */
		list = _protocolFile.buildSetMarkUp("species");
		speciesDic = new ArrayList<String>(list.size());
		for (Element aChild : list)
			speciesDic.add(aChild.getAttributeValue("name"));

		/* Build the dictionary of solvers ___________________________ */
		list = _protocolFile.buildSetMarkUp("solver");
		solverDic = new ArrayList<String>(list.size());
		solverList = new DiffusionSolver[list.size()];
		for (Element aChild : list)
			solverDic.add(aChild.getAttributeValue("name"));

	}

	/**
	 * Create the world properties e.g. system size and boundary conditions.
	 * Currently you can create zero-flow, cyclic, constant and dilution
	 * boundary ; to create on on your own, create a new class extending the
	 * abstract class "BoundaryCondition"
	 */
	public void createWorld() {
		// Creation of the world
		try {
			System.out.print("\t World: \n");
			XMLParser parser = new XMLParser(_protocolFile.getChildElement("world"));
			world = new World();
			world.init(this, parser);

			// now set the bulk concentrations if it is needed
			if (useBulkFile) recreateBulkConditions();

			System.out.println("\t done");
		} catch (Exception e) {
			LogFile.writeLog("Simulator.createWorld() : "+e);
			System.exit(-1);
		}
	}

	/**
	 * Get the bulk concentrations from the input file and assign them to the
	 * current bulks
	 * 
	 * @added 23.1.2009
	 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
	 */
	public void recreateBulkConditions() throws Exception {
		String bulkName;
		int soluteIndex;
		String soluteName;
		XMLParser simulationRoot = new XMLParser(bulkFile.getChildElement("simulation"));

		for (Element aBulkMarkUp : simulationRoot.buildSetMarkUp("bulk")) {

			XMLParser aBulkRoot = new XMLParser(aBulkMarkUp);
			bulkName = aBulkRoot.getAttributeStr("name");

			// check to make sure the bulk exists
			if (!world.containsBulk(bulkName))
				throw new Exception("Bulk "+bulkName+" is not specified in protocol file");

			LogFile.writeLog("\t\tInitializing bulk '"+bulkName+"' from input file.");

			Bulk thisBulk = world.getBulk(bulkName);

			// now set the solutes within this bulk
			for (Element aSoluteMarkUp : aBulkRoot.buildSetMarkUp("solute")) {
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
	 * Create all soluble species as spatial grids. The protocol File specifies
	 * the solute name, its diffusivity, eventually the connected bulks and
	 * their concentrations
	 */
	public void createSolutes() {
		System.out.print("\t Solutes: \n");
		try {
			int iSolute = 0;
			for (Element aSoluteMarkUp : _protocolFile.buildSetMarkUp("solute")) {
				XMLParser aSoluteRoot = new XMLParser(aSoluteMarkUp);
				soluteList[iSolute] = new SoluteGrid(this, aSoluteRoot);
				LogFile.writeLog("\t\t"+soluteList[iSolute].getName()+" ("+soluteList[iSolute].soluteIndex+")");
				iSolute++;
			}
			System.out.println("\t done");
		} catch (Exception e) {
			LogFile.writeLog("Simulator.createSolutes() : error met " + e);

			System.exit(-1);
		}
	}

	/**
	 * Create all reactions described in the protocol file
	 */
	public void createReactions() throws Exception {
		Reaction aReaction;

		LogFile.writeLog("\t Reactions: \n");

		int iReaction = 0;
		for (Element aReactionMarkUp : _protocolFile.buildSetMarkUp("reaction")) {
			XMLParser aReactionRoot = new XMLParser(aReactionMarkUp);
			aReaction = (Reaction) aReactionRoot.instanceCreator("simulator.reaction");
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
	 * Create the solvers for the diffusion-reaction-system described in the
	 * protocol file.
	 */
	public void createSolvers() throws Exception {
		System.out.print("\t Solvers: \n");
		try{
		for (Element aSolverMarkUp : _protocolFile.buildSetMarkUp("solver")) {
			// Initialise the XML parser
			XMLParser parser = new XMLParser(aSolverMarkUp);

			// Create the solver,initialise it and register it
			DiffusionSolver aSolver =
				(DiffusionSolver) parser.instanceCreator("simulator.diffusionSolver");
			aSolver.init(this, parser);
			aSolver.register();
			System.out.println("\t\t"+aSolver.solverName+" ("+aSolver.solverIndex+")");
		}

		System.out.println("\t done");
		}catch(Exception e){LogFile.writeLog("Error in Simulator.creatSolvers(): "+e);}
	}

	/**
	 * Create the agentGrid, the species and the agents
	 */
	public void createSpecies() throws Exception {
		//try {

			try { // First stage
				// Create the Species (and the progenitor) and register it
				System.out.print("\t Species: \n");
				for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) {
					Species aSpecies = new Species(this, new XMLParser(aSpeciesMarkUp)); 
					speciesList.add(aSpecies);
					LogFile.writeLog("\t\t"+aSpecies.speciesName+" ("+aSpecies.speciesIndex+")");
				}
				System.out.print("\t done\n");
			} catch (Exception e) {
				LogFile.writeLog("Error in Simulator.createSpecies() first stage");
			}

			try { // Second stage
				// Create the agent grid
				System.out.print("\t Agent Grid: \n");
				XMLParser parser = new XMLParser(_protocolFile.getChildElement("agentGrid"));
				agentGrid = new AgentContainer(this, parser, agentTimeStep);
				is3D = agentGrid.is3D;
				System.out.print("\t done\n");


				// Finalise the initialisation of the progenitor
				System.out.print("\t Species progenitor: \n");
				for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) {
					parser = new XMLParser(aSpeciesMarkUp);
					getSpecies(parser.getAttribute("name")).getProgenitor().initFromProtocolFile(this,
							parser);
					//sonia: creating a list with the plasmid names which will be used afterwards to write the agentSum report

					if(multiEpi){
						if(parser.getAttribute("class").equals("MultiEpisome")){
							String plName = parser.getAttribute("name");
							Double scanSpeed = parser.getParamDbl("scanSpeed");
							plasmidList.add(plName);
							scanSpeedList.add(scanSpeed);

						}
					}
				}
				System.out.print("\t done\n");
			} catch (Exception e) {
				LogFile.writeLog("Error in Simulator.createSpecies() second stage");
			}
			
			//try { // Third stage
				// Create the population
				System.out.print("\t Species populations: \n");

				if (useAgentFile){
					recreateSpecies();
					
					//sonia: 20-07-09 
					//I've added the line of code below, so that it is possible to create new agents from the protocol file when restarting
					//from a previous simulation (environment) using a new protocol file
					checkAgentBirth();

				}else{
					checkAgentBirth();
				}
				System.out.print("\t done\n");
			//} catch (Exception e) {
				//LogFile.writeLog("Error in Simulator.createSpecies() third stage");
			//}
		//} catch (Exception e) {
		//	LogFile.writeLog("Simulator.createSpecies(): error met: "+e);}
	}

	/**
	 * Create the agents from a file describing individually all the agents
	 * species by species
	 * @throws Exception
	 * @author ad
	 */
	public void recreateSpecies() throws Exception {
		int spIndex, counterSpecies;
		SpecialisedAgent progenitor;
		XMLParser simulationRoot = new XMLParser(agentFile.getChildElement("simulation"));
		counterSpecies = 0;
		
		for (Element aSpeciesMarkUp : simulationRoot.buildSetMarkUp("species")) {

			XMLParser aSpeciesRoot = new XMLParser(aSpeciesMarkUp);
			spIndex = getSpeciesIndex(aSpeciesRoot.getAttribute("name"));
			
			// check consistency with protocol file declarations
			boolean isConsistent = (speciesList.get(counterSpecies).speciesName.equals(aSpeciesRoot
					.getAttributeStr("name")));
			if (!isConsistent) throw new Exception(
			"Agent input file is inconsistent with protocol file: ");
			
			// Process agents description
			String dataSource = aSpeciesRoot.getElement().getContent(0).toString();
			String[] allAgentData = null;
			try {
				allAgentData = dataSource.split(";\n");
				// this removes the '[Text: \n' line from the first item in the list
				// so that all agents will be treated equally
			} catch (Exception e){
				LogFile.writeLog("Simulator.recreateSpecies() : problem splitting up data");}
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
	 * Test if new bacteria should be created in the system
	 */
	public void checkAgentBirth() {
		XMLParser parser;
		int spIndex;
		boolean creatingAgents = false;

		for (Element aSpeciesMarkUp : _protocolFile.buildSetMarkUp("species")) {
			parser = new XMLParser(aSpeciesMarkUp);
			spIndex = getSpeciesIndex(parser.getAttribute("name"));

			for (Element aInitMarkUp : parser.buildSetMarkUp("initArea")) {
				parser = new XMLParser(aInitMarkUp);
				if (SimTimer.isDuringNextStep(parser.getParamTime("birthday"))) {
					this.agentGrid.agentIter = agentGrid.agentList.listIterator();
					speciesList.get(spIndex).createPop(parser);
					creatingAgents = true;
				}
			}
		}
		if (creatingAgents) agentGrid.relaxGrid();
	}

	/* _____________________________________________________________________ */

	/**
	 * Generate a report with concentrations of solutes and biomass on the
	 * default grid and a report with the exhaustive description of all agents
	 * Each report is a new file with a new index
	 */
	public void writeReport() {
		// Update saving counters and file index
		_lastOutput = SimTimer.getCurrentTime();
		int currentIter = SimTimer.getCurrentIter(); // bvm added 26.1.2009

		// first restart log file to avoid non-write trouble
		LogFile.reopenFile();

		try {
			/* Grids and environment ______________________________________ */
			// env_State
			result[0].openFile(currentIter);
			// env_Sum
			result[1].openFile(currentIter);
			// bvm added 16.12.08


			//sonia:chemostat
			if(Simulator.isChemostat){
				//sonia:chemostat
				//I've modified refreshBiofilmGrids()
				soluteList[0].getDomain().refreshBioFilmGrids();
			}else{


				// output the biofilm thickness data
				//sonia 12.10.09				

				double [] intvals;
				StringBuffer value = new StringBuffer();
				for (Domain aDomain : world.domainList) {

					//double _resolution = aDomain._resolution;
					aDomain.refreshBioFilmGrids();
					intvals = aDomain.getInterface();

					value.append("<thickness domain=\""+aDomain.domainName+"\" unit=\"um\">\n");
					value.append("\t<mean>"+(ExtraMath.mean(intvals))+"</mean>\n");
					value.append("\t<stddev>"+(ExtraMath.stddev(intvals))+"</stddev>\n");
					value.append("\t<max>"+(ExtraMath.max(intvals))+"</max>\n");
					value.append("</thickness>\n");

				}

				result[0].write(value.toString());
				result[1].write(value.toString());

			}

			// Add description of each solute grid
			for (SoluteGrid aSG : soluteList) {
				aSG.writeReport(result[0], result[1]);
			}

			// Add description of each reaction grid
			for (Reaction aReac : reactionList) {
				aReac.writeReport(result[0], result[1]);
			}

			// Add description of each species grid
			agentGrid.writeGrids(this, result[0], result[1]);

			// Add description of total biomass
			for (Domain aDomain : world.domainList) {
				aDomain.refreshBioFilmGrids();
				aDomain.getBiomass().writeReport(result[0], result[1]);
			}

			// Add description of bulks
			for (Bulk aBulk : world.bulkList) {
				aBulk.writeReport(result[1]);
			}

			result[0].closeFile();
			result[1].closeFile();

		} catch (Exception e) {
			LogFile.writeError("System description of grids failed:"+e.getMessage(),
			"Simulator.writeReport()");
		}

		try {
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

		} catch (Exception e) {
			LogFile.writeError("System description of agents failed:"+e.getMessage(),
			"Simulator.writeReport()");
		}

	}


	/* ____________________ ACCESSORS & MUTATORS ___________________________ */
	/**
	 * Find a Species on the basis of its nickname
	 * 
	 * @param aSpeciesName
	 * @return the speciesIndex
	 */
	public int getSpeciesIndex(String aSpeciesName) {
		return speciesDic.indexOf(aSpeciesName);
	}

	public Species getSpecies(String aSpeciesName) {
		return speciesList.get(getSpeciesIndex(aSpeciesName));
	}

	public int getSoluteIndex(String aSoluteName) {
		return soluteDic.indexOf(aSoluteName);
	}

	public SoluteGrid getSolute(String aSoluteName) {
		return soluteList[getSoluteIndex(aSoluteName)];
	}

	public int getReactionIndex(String aReactionName) {
		return reactionDic.indexOf(aReactionName);
	}

	public Reaction getReaction(String aReactionName) {
		return reactionList[getReactionIndex(aReactionName)];
	}

	public int getSolverIndex(String aSolverName) {
		return solverDic.indexOf(aSolverName);
	}

	public DiffusionSolver getSolver(String aSolverName) {
		int solInd = getSolverIndex(aSolverName);
		if (solInd >= 0)
			return solverList[solInd];
		else
			return null;
	}

	public int getParticleIndex(String particleName) {
		return particleDic.indexOf(particleName);
	}

	/* _____________________________ GUI ________________________________ */

	public void play() {
	}

	public void pause() {
	}

	public void stop() {
	}

	/**
	 * Check whether we will initialize from agent and bulk input files
	 * If the "restartPreviousRun" param is true, this will also set
	 * the correct files for reading in the last state
	 * 
	 * @param protocolFile
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

	public String getResultPath() {
		return _resultPath;
	}

}
