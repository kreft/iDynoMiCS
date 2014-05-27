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

import java.util.*;

import simulator.agent.*;
import simulator.agent.zoo.MultiEpiBac;

import simulator.detachment.*;
import simulator.diffusionSolver.Solver_pressure;

import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.reaction.Reaction;
import simulator.SpatialGrid;

import utils.ResultFile;
import utils.XMLParser;
import utils.LogFile;
import utils.ExtraMath;

/**
 * \brief Class to store all the agents, call them, and manage shoving/erosion of located agents
 * 
 * Class to store all the agents, call them, and manage shoving/erosion of located agents
 * 
 * @author Andreas Doetsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sonia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class AgentContainer 
{

	/**
	 * Computational domain to which this grid is assigned
	 */ 
	public Domain domain;
	
	/**
	 * Local copy of the simulation object used to create the conditions specified in the protocol file
	 */
	public Simulator mySim;

	/**
	 * Container for all agents (even the non located ones)
	 */
	public LinkedList<SpecialisedAgent> agentList;
	
	/**
	 * Iterator for all agents in this grid
	 */
	public ListIterator<SpecialisedAgent> agentIter;

	/**
	 * Temporary containers used to store agents who will be added or removed.  Visibility public so that it can be accessed from LocatedGroup in killAll()
	 */
	public LinkedList<SpecialisedAgent> _agentToKill = new LinkedList<SpecialisedAgent>();

	/**
	 * Array of SpatialGrids - one for each species in the simulation
	 */ 
	private SpatialGrid[] _speciesGrid;
	
	/**
	 * Number of grid elements in the x direction
	 */
	private int _nI;
	
	/**
	 * Number of grid elements in the y direction
	 */
	private int _nJ;
	
	/**
	 * Number of grid elements in the z direction
	 */
	private int _nK;
	
	/**
	 * Total number of grid elements in this container
	 */
	private int _nTotal;
	
	/**
	 * Resolution of the grid. Specified in the XML protocol file
	 */
	private double _res;
	
	
	/**
	 *	Grid to hold the agents in this container. This grid holds groups of agents, within a LocatedGroup object 
	 */
	private LocatedGroup[] _grid;
	
	/**
	 * 3D array that captures erosion in the agent grid
	 */
	protected double[][][] _erosionGrid;
	
	/**
	 * Boolean noting whether the grid is 3D or 2D
	 */ 
	public boolean is3D;
	
	/**
	 * Sets the allowed cut off for moving agents
	 */
	private final double SHOV_FRACTION;

	/**
	 * Allows you to define a smaller timestep for agent behaviors and interactions. Read in by Simulator class. This is a local copy 
	 */
	public final double AGENTTIMESTEP;
	
	/**
	 * Sets the maximal number of iterations allowed to reach the final agent positions
	 */
	private final int MAXITER;
	
	/**
	 * Whether shoving motion is applied to both agents or only to one when two agents overlap
	 */
	private final boolean MUTUAL;
	
	/**
	 * Boolean to select erosion method: true is shrinkOnBorder, false is removeOnBorder
	 */
	private final boolean EROSIONMETHOD;
	
	/**
	 * Whether biofilm pieces that are no longer connected to the substratum are removed
	 */
	private final boolean DOSLOUGHING;
	
	/**
	 * Limit to agent number that will cause the simulation to cease
	 */
	private final int maxPopLimit;
	
	/**
	 * Grid used to store pressure, created if specified in the protocol file
	 */
	private SoluteGrid _pressure;
	
	/**
	 *  Tally of mass to be removed in removeOnBorder, or cells in chemostat dilution (see agentsFlushedAway)
	 */
	double tallyVariable = 0; 
	
	/**
	 * Number of shoving iterations performed
	 */
	int shovIter;
	
	/**
	 * Limit on number of cells that will be moved in shoving
	 */
	int shovLimit; 
	
	/**
	 * Maximum number of shove iterations that can be performed in a step
	 */
	int maxShoveIter; 
	
	/**
	 * Number of agents that have been moved due to shoving
	 */
	int nMoved;
	
	/**
	 * Number of cells that perform a 'relative' move under shoving
	 */
	double tMoved; 
	
	/**
	 * Calculation of a move that is to be performed due to influence of pressure 
	 */
	double deltaMove;

	/**
	 * Solver used for modelling detachment
	 */
	public LevelSet _levelset;

	/**
	 * Calculated factor noting the influence of dilution on the agent grid. Sonia Martins 23.11.09
	 */
	private double Dfactor;

	/**
	 * \brief Creates the agent grid in which all agents in the biofilm simulations are stored, as well as erosion and species grids
	 * 
	 * All agents are stored within a grid, called the agentGrid. This is similar to that for solutes but serves a purpose unique 
	 * to agents. This constructor initialises this grid with parameters taken from the XML file
	 * 
	 * @param aSimulator	The simulation object used to simulate the conditions specified in the protocol file
	 * @param root	The agentGrid markup from the XML file
	 * @param agentTimeStep	The agent time step as specified in the XML protocol file
	 * @throws Exception	Exception thrown should this data not be present in the protocol file
	 */
	public AgentContainer(Simulator aSimulator, XMLParser root, double agentTimeStep) throws Exception 
	{
		// Read FINAL fields
		SHOV_FRACTION = root.getParamDbl("shovingFraction");
		MAXITER = root.getParamInt("shovingMaxIter");
		MUTUAL = root.getParamBool("shovingMutual");

		
		if (root.getParam("erosionMethod") == null)
			EROSIONMETHOD = true;
		else
			EROSIONMETHOD = root.getParamBool("erosionMethod");

		if (root.getParam("sloughDetachedBiomass") == null)
			DOSLOUGHING = true; // default to carry out sloughing
		else
			DOSLOUGHING = root.getParamBool("sloughDetachedBiomass");
			
		if (root.getParamBool("maxPopLimit") == null)
			maxPopLimit = 0;
		else
			maxPopLimit = root.getParamInt("maxPopLimit");

		// Now deal with the agent timestep
		if (Double.isNaN(agentTimeStep))
		{
			AGENTTIMESTEP = SimTimer.getCurrentTimeStep();
			LogFile.writeLog("Using global timestep of "+AGENTTIMESTEP+" for agentTimeStep");
		}
		else
		{
			AGENTTIMESTEP = agentTimeStep;
			if (AGENTTIMESTEP > SimTimer.getCurrentTimeStep()) 
			{
				LogFile.writeLog("ERROR: agentTimeStep in agentGrid markup MUST be "+
						"less than or equal to the global timestep\n"+
						"\tagentTimeStep was given as: "+AGENTTIMESTEP+"\n"+
						"\tglobal time step is currently: "+SimTimer.getCurrentTimeStep());
				throw new Exception("agentTimeStep too large");
			}
			LogFile.writeLog("Agent time step is... " + agentTimeStep);
		}

		// Now set the domain where this container is defined
		domain = aSimulator.world.getDomain(root.getParam("computationDomain"));
		mySim = aSimulator;
		
		agentList = new LinkedList<SpecialisedAgent>();
		agentIter = agentList.listIterator();
		// Optimised the resolution of the grid used to sort located agents
		checkGridSize(aSimulator, root);

		// Now initialise the padded grid
		createShovGrid(aSimulator);
		
		// Initialise spatial grid used to display species distribution
		createOutputGrid(aSimulator);

		if (Simulator.isChemostat)
		{
			LogFile.writeLog("Chemostat volume is "+ExtraMath.cube(_res)+" micrometers cubed");
		}
		else
		{
			// initialise the pressure grid, if there is one
			if (aSimulator.getSolver("pressure")!=null)
				_pressure = ((Solver_pressure) aSimulator.getSolver("pressure")).getPressureGrid();

			// Initialise the levelset solver used for detachment
			_levelset = LevelSet.staticBuilder(root.getChild("detachment"), this);

			LogFile.writeLog(" " + _nTotal + " grid elements, resolution: " + _res
					+ " micrometers");
		}


	}

	/* ___________________ STEPPERS ________________________________ */

	/**
	 * \brief Iterates through each grid cell on the agent grid, stepping all agents that are within that grid space 
	 * 
	 * Iterates through each grid cell on the agent grid, stepping all agents that are within that grid space
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 */
	public void step(Simulator aSim) {

		SpecialisedAgent anAgent;
		int nDead, nAgent, nBirth;

		/* STEP AGENTS ________________________________________________ */
		LogFile.chronoMessageIn();
		Collections.shuffle(agentList, ExtraMath.random);

		// record values at the beginning
		nDead = 0;
		nBirth = 0;
		nAgent = agentList.size();
		double dt = 0;
		double elapsedTime = 0;
		double globalTimeStep = SimTimer.getCurrentTimeStep();
		// for the local time step, choose the value according to which is best
		double localdt = Math.min(AGENTTIMESTEP,globalTimeStep);

		double nAgent0 = agentList.size();
		// Apply a shorter time step when visiting all the agents

		while (elapsedTime < globalTimeStep) {
			// by default use the saved agent timestep
			dt = localdt;


			// check for a smaller dt (usually for the last iterate)
			if (dt > (globalTimeStep-elapsedTime))
				dt = globalTimeStep-elapsedTime;

			elapsedTime += dt;		

			/* Step all the agents */
			SimTimer.setCurrentTimeStep(dt);

			//sonia:chemostat
			if(Simulator.isChemostat){
				//sonia: bypass bacteria movement according to pressure field
			}else{
				followPressure();
			}


			for (agentIter = agentList.listIterator(); agentIter.hasNext();) {
				anAgent = agentIter.next();
				anAgent.step();
			}

			Collections.shuffle(agentList, ExtraMath.random);

			if (Simulator.isChemostat){
				agentFlushedAway(dt);		
			}


			// Add and remove agents
			nBirth += agentList.size() - nAgent;

			//sonia 26.04.2010
			//commented out removeAllDead
			// this call is now made at the end of the step in Simulator
			//nDead += removeAllDead();
			nAgent = agentList.size();

			// NOW DEAL WITH DEATH IN THIS AGENT TIMESTEP
			// REMOVE THESE FROM THE GRID IF DEAD
			// MUST BE DONE SO THAT THESE DO NOT AFFECT SHOVING
			for(SpecialisedAgent aDeathAgent: _agentToKill)
			{
				if (aDeathAgent.isDead) 
				{
					//nDead++;
					// KA - removed the count here, as the count of dead cells was wrong - we were recounting these with every
					// agent timestep. We should only be counting them at the simulation timestep
					// However they need to remain in the _agentToKill list until this is emptied at the correct output period
					agentList.remove(aDeathAgent);
					removeLocated(aDeathAgent);
				}
			}

			// Apply moderate overlap relaxation

			//sonia:chemostat
			//11-06-09
			if(!Simulator.isChemostat){
				shoveAllLocated(false, true, 15, 1, 1);
			}
			
		}

		SimTimer.setCurrentTimeStep(globalTimeStep);

		//sonia 11.10.2010 - implementing HGT step after agents' division and shoving

		/*	double elapsedHGTtime=0;
		double hgtTimeStep =0;
		double maxScan =0;
		double currentTime = SimTimer.getCurrentTime();

		if (Simulator.multiEpi){

			elapsedHGTtime = currentTime;

			maxScan = mySim.scanSpeedList.get(0);
			for (int i=0; i< mySim.scanSpeedList.size(); i++ ){
				if (mySim.scanSpeedList.get(i)> maxScan){
					maxScan = mySim.scanSpeedList.get(i);
				}
			}

			hgtTimeStep = 1/ maxScan;
			System.out.println("hgt time step is " + hgtTimeStep);

			while (elapsedHGTtime < (globalTimeStep + currentTime)) {
				//	System.out.println("elapsedHGTtime is " + elapsedHGTtime);

				for (agentIter = agentList.listIterator(); agentIter.hasNext();) {
					anAgent = agentIter.next();
					anAgent.HGTstep(elapsedHGTtime);
				}

				elapsedHGTtime += hgtTimeStep;

			}
		}*/

		// KA - MOVED OUTPUT OF AGENTS STEPPED / DEAD / BORN FROM HERE TO LATER, SO THAT WE CAN INCLUDE ERODED CELLS IN THE
		// COUNT OF DEAD CELLS
				
		/* MECHANICAL INTERACTIONS _____________________________________ */

		//sonia:chemostat
		//11-06-09
		if(Simulator.isChemostat){
			// bypass the shoving, detachment and erosion processes 

		}else{

			//sonia 26.04.2010
			//care as been take so that the death agents are removed from the 
			// _agentList preventing their participation in the shoving process 			
			// spring and then shove only particles
			shoveAllLocated(false, true, MAXITER, 1, 1);
			LogFile.writeLog(nMoved + "/" + agentList.size() + " after " + shovIter
					+ " shove iterations");



			// EROSION & DETACHMENT _________________________________________ */
			// Refresh the space occupation map (-1:outside, 0:carrier,1:biofilm, 2:liquid, 3:bulk)
			refreshGroupStatus();

			// Rebuild the border of the biofilm and compute erosion-time for the
			// whole biofilm

			_levelset.refreshBorder(true, mySim);
			_levelset.computeLevelSet(mySim);

			// On grid elements on the border apply a probabilistic erosion
			// Rob Feb 2011: added borderErosion(), an alternative to erodeBorder(), which
			// removes whole agents in a deterministic way
			if (EROSIONMETHOD)
				shrinkOnBorder();
			else
				// KA - Check added after self-attachment, as it is possible that if the timestep is much less than the input rate, 
				// the grid may contain no cells for the first few steps
				if(this.agentIter.hasNext())
				{
					try {
						removeOnBorder(this);
					} catch (Exception e) {
						System.out.println("At AgentContainer:removeOnBorder error met :" + e);
						e.printStackTrace();
						System.exit(-1);
					}
					// mark biomass connected to the carrier and remove any non-connected portions
					if (DOSLOUGHING) {
						refreshGroupStatus();
						markForSloughing();
					}
				}

			LogFile.chronoMessageOut("Detachment");
			
				
		}
		
		// OUTPUT THE COUNT STATISTICS
		LogFile.chronoMessageOut("Agents stepped/dead/born: " + nAgent0 + "/"
				+ _agentToKill.size() + "/" + nBirth);

		
		nAgent = agentList.size();
		if (maxPopLimit > 0 && nAgent >= maxPopLimit)
			aSim.continueRunning = false;
	}


	/**
	 * \brief Compute pressure field and apply resulting advection movement to affected agents
	 * 
	 * Compute pressure field and apply resulting advection movement to affected agents
	 */
	public double followPressure() 
	{
		double moveMax;

		// Find a solver for pressure field and use it
		if (mySim.getSolver("pressure")==null) return 0;

		// don't use the pressure if it's not active
		if (!mySim.getSolver("pressure").isActive()) return 0;

		LogFile.writeLog("Doing pressure calculations.");


		// get local timestep (which was set in the step() routine calling this one)
		double dt = SimTimer.getCurrentTimeStep();


		// Solve for pressure field
		mySim.getSolver("pressure").initAndSolve();
		_pressure = ((Solver_pressure) mySim.getSolver("pressure")).getPressureGrid();


		// copy calculated pressure field to the solute list
		// (allows easy output of pressure field)
		mySim.getSolute("pressure").setGrid(_pressure.getGrid());

		// Determine local advection speed
		moveMax = 0;
		for (LocatedGroup aGroup : _grid)
			moveMax = Math.max(moveMax, aGroup.computeMove(_pressure,AGENTTIMESTEP));

		// bvm 04.03.09: new method to address any high velocities:
		// use smaller local timesteps to keep the movement under control
		double dtlocal = dt;
		int itlocal = 1;
		while (dtlocal > this._res/moveMax) {
			// if the move takes an agent farther than one grid element,
			// apply scaling factor until move is within limit
			dtlocal /= 10.;
			itlocal *= 10;
		}
		if (itlocal > 1) {
			LogFile.writeLog("PRESSURE MOVEMENT HAS LOCAL TIMESTEP "
					+dtlocal+" ("+itlocal+" iterations)");
		}

		// scale movement vectors based on new, smaller timestep and apply
		// the movement to each agent in each group
		double alpha = dtlocal/dt;
		for (LocatedGroup aGroup : _grid)
			aGroup.addMoveToAgents(alpha);

		// now apply the scaled agent movements to each agent
		deltaMove = 0;
		for (int i=0; i<itlocal; ++i) {
			agentIter = agentList.listIterator();
			while (agentIter.hasNext())
				deltaMove += agentIter.next().move();
		}

		return deltaMove;
	}	


	/**
	 * \brief Solve spatial spreading through application of shoving (acts only on located agents)
	 * 
	 * Solve spatial spreading through application of shoving (acts only on located agents)
	 * 
	 * @param fullRelax	Boolean noting whether a full relax of the grid is applied
	 * @param shoveOnly	Boolean noting whether to apply shoving only or to also apply agent pushing
	 * @param maxShoveIter	The maximum number of shoving iterations that should be applied to find a new position
	 * @param gainMin	Specified minimum amount of gain to be found
	 * @param gainMax	Specified maximum amount of gain to be found
	 */
	public void shoveAllLocated(boolean fullRelax, boolean shoveOnly,
			double maxShoveIter, double gainMin, double gainMax) {

		if (fullRelax)
			maxShoveIter = MAXITER * 5;

		nMoved = shovLimit = Math.max(1,
				(int) (agentList.size() * SHOV_FRACTION));
		shovIter = 0;
		while ((nMoved >= shovLimit) & (shovIter++ < maxShoveIter)) {
			performMove(shoveOnly, false, 1);
		}
	}

	/**
	 * \brief Used during initialisation to start from a coherent state
	 * 
	 * Used during initialisation to start from a coherent state
	 */
	public void relaxGrid() 
	{

		//sonia:chemostat
		if(!Simulator.isChemostat)
		{
			Collections.shuffle(agentList, ExtraMath.random);
			shoveAllLocated(true, true, MAXITER / 2, 0.1, 0.25);
			shoveAllLocated(true, true, MAXITER / 2, 0.1, 1);
		}
	}


	/**
	 * \brief Moves an agent as a result of shoving or pushing due to growth
	 *
	 * Moves an agent as a result of shoving or pushing due to growth
	 * 
	 * @param pushOnly	Boolean noting whether to apply shoving only or to also apply agent pushing
	 * @param isSynchro
	 * @param gain
	 */
	protected double performMove(boolean pushOnly, boolean isSynchro,
			double gain) {
		SpecialisedAgent anAgent;
		nMoved = 0;
		tMoved = 0;
		double nMoved2 = 0;

		for (agentIter = agentList.listIterator(); agentIter.hasNext();) {
			// Compute movement, deltaMove is relative movement
			anAgent = agentIter.next();
			deltaMove = anAgent.interact(MUTUAL, pushOnly, !isSynchro, gain);

			tMoved += deltaMove;
			nMoved += (deltaMove >= 0.1 * gain ? 1 : 0);
			nMoved2 += (deltaMove >= 0.1 ? 1 : 0);
		}


		if (!isSynchro)
			return nMoved2;

		for (agentIter = agentList.listIterator(); agentIter.hasNext();) {
			// Compute movement, deltaMove is relative movement
			anAgent = agentIter.next();
			deltaMove = anAgent.move();

			tMoved += deltaMove;
			nMoved += (deltaMove >= 0.1 * gain ? 1 : 0);
			nMoved2 += (deltaMove >= 0.1 ? 1 : 0);

			if (anAgent.isDead){
				anAgent.death = "invalidMove";
				//_agentToKill.add(anAgent);
				//sonia 26.04.2010
				//added agentList.remove(anAgent);
				//agentList.remove(anAgent);
			}
		}
		//sonia 26.04.2010
		//commented out removeAllDead()
		//removeAllDead();
		return nMoved2;
	}


	/**
	 * \brief Refresh the space occupation map as agents may have moved
	 * 
	 * Refresh the space occupation map as agents may have moved. Each grid square has a status: -1:outside, 0:carrier,1:biofilm, 
	 * 2:liquid, 3:bulk
	 */
	protected void refreshGroupStatus() {
		for (int index = 0; index < _nTotal; index++)
			_grid[index].refreshElement();
	}

	/**
	 * \brief Explore grid cells around the current one and returns all agents that may shove this agent.
	 * 
	 * Explore grid cells around the current one and returns all agents that may shove this agent. Includes grid cells on the other 
	 * side of the cyclic boundary
	 * 
	 * @param index	The integer index of the grid square on the agent grid
	 * @param range	maximal range to screen for shoving agents
	 * @param nbList: the list of located agents
	 */
	public void getPotentialShovers(int index, double range,
			LinkedList<LocatedAgent> nbList) {
		LocatedGroup aGroup;
		int radius = Math.max(1, (int) Math.floor(range / this._res));

		nbList.clear();

		for (int i = -radius; i <= radius; i++) {
			if (_grid[index].moveX(i) == null)
				continue;
			for (int j = -radius; j <= radius; j++) {

				if (!is3D) {
					aGroup = _grid[index].moveX(i).moveY(j);

					if (aGroup != null)
						nbList.addAll(aGroup.group);

				} else {
					if (_grid[index].moveX(i).moveY(j) == null)
						continue;
					for (int k = -radius; k <= radius; k++) {
						aGroup = _grid[index].moveX(i).moveY(j).moveZ(k);
						if (aGroup != null)
							nbList.addAll(aGroup.group);
					}
				}
			}
		}
	}

	/* ________________ TOOLS:GRID, MAP & TREE MANAGEMENT __________________ */

	/**
	 * \brief Registers the birth of an agent and adds this to the agent grid
	 * 
	 * Registers the birth of an agent and adds this to the agent grid
	 * 
	 * @param anAgent	New agent to add to the agent grid
	 */
	public void registerBirth(SpecialisedAgent anAgent) 
	{
		// Add the agent to agentList
		agentIter.add(anAgent);

		// Add the agent on the grid
		if (anAgent instanceof LocatedAgent) {
			LocatedAgent aLoc = (LocatedAgent) anAgent;
			try {
				if(Simulator.isChemostat){
					_grid[0].add(aLoc);
				}else{
					int index = getIndexedPosition(aLoc.getLocation());
					if (!Double.isNaN(index))
						_grid[index].add(aLoc);
				}
			} catch (Exception e) {
				LogFile.writeLog("Error:Failed to add an agent on the grid");
			}
		}

	}

	/**
	 * \brief Register the death of an agent, by adding it to a list of agents that will be removed from the simulation
	 * 
	 * Register the death of an agent, by adding it to a list of agents that will be removed from the simulation
	 * 
	 * @param anAgent	SpecialisedAgent object that will be removed from the simulation
	 */
	public void registerDeath(SpecialisedAgent anAgent) 
	{
		if(!_agentToKill.contains(anAgent))
			_agentToKill.add(anAgent);

	}

	/**
	 * \brief Iterates through the _agentToKill list and removes these agents from the grid. Returns the number of dead cells removed
	 * 
	 * Iterates through the _agentToKill list and removes these agents from the grid. These cells should be assumed to be dead
	 * 
	 * @return	Number of dead cells removed from the simulation
	 */
	public int removeAllDead() {
		int nDead = 0;

		ListIterator<SpecialisedAgent> iter = _agentToKill.listIterator();
		SpecialisedAgent anAgent;

		while (iter.hasNext()) {
			anAgent = iter.next();
			if (anAgent.isDead) {
				nDead++;
				iter.remove();
				agentList.remove(anAgent);
				removeLocated(anAgent);
			}
		}

		_agentToKill.clear();
		return nDead;
	}

	/**
	 * \brief Removes a number of agents from the system according to the dilution set for the chemostat.
	 * 
	 * Removes a number of agents from the system according to the dilution set for the chemostat. The global time step should be set 
	 * to 0.10*(1/D) so that around 10% of the agents will be removed from the system in each iteration. Remember, agents stand for all 
	 * type of particles that can be removed (aka deleted) from the system, from bacteria to eps.
	 *
	 * @param agentTimeStep - this should be the same as the global timeStep or lower
	 */
	public void agentFlushedAway(double agentTimeStep){

		// after having shuffled the list (during the step()) with all the agents we are now ready to kill
		// agents according to the dilution value read from the Bulk class

		//sonia:2.03.2010 Just in case, let's do it again
		Collections.shuffle(agentList, ExtraMath.random);

		//double randNum;

		for (AllBC aBC : domain.getAllBoundaries()){
			if (aBC.hasBulk()){
				Bulk aBulk = aBC.getBulk();
				if(aBulk.getName().equals("chemostat")){
					Dfactor = aBulk._D;
				}
			}
		}
		
		// Rob 9/6/11: simplified and removed any effects which might be caused by rounding
		//agentsToDilute =  Math.round(Dfactor*(agentTimeStep)* agentList.size());
		//pDying = agentsToDilute / agentList.size();
//		pDying = Dfactor*agentTimeStep;
//
//		for ( int i=0; i< agentList.size(); i++){
//			randNum = ExtraMath.getUniRand();
//			if(randNum < pDying){
//				agentList.get(i).isDead = true;
//				agentList.get(i).death = "dilution";
//				agentList.get(i).die(false);
//				dead ++;
//			}
//		}
		int agentsToDilute = 0;
		if (EROSIONMETHOD){
//			agentsToDilute = (int) Math.round(Dfactor*agentTimeStep*agentList.size());
			// Rob 08/02/2012: Instead of simply rounding, we now use the non-integer part in a
			// stochastic way. This means that washout of populations is possible when using a
			// rate and an agent time-step whose product is less than 0.5
//			double temp = Dfactor*agentTimeStep*agentList.size();
//			if (Math.floor(temp) + Math.random() < temp)
//				agentsToDilute = (int) Math.floor(temp);
//			else
//				agentsToDilute = (int) Math.ceil(temp);			
			// Rob 13/09/2012: this method is now changed t@param aSim	The simulation object used to simulate the conditions specified in the protocol fileo include a running tally of the non-integer
			// remainder (as in removeOnBorder in biofilms)
			double temp = Dfactor*agentTimeStep*agentList.size() + tallyVariable;
			agentsToDilute = (int) Math.floor(temp);
			tallyVariable = temp%1;
			
		} else{
			// Rob 28/11/2011: Added this so we can keep the population at or below 1000
			// EROSIONMETHOD is normally used to decide between biofilm functions so this shouldn't be a problem
			agentsToDilute = Math.max(agentList.size() - 1000, 0);
		}
		
		for (int i=0; i<agentsToDilute; i++){
			agentList.get(i).isDead = true;
			agentList.get(i).death = "dilution";
			agentList.get(i).die(false);
		}

		SpecialisedAgent anAgent;
		agentIter = agentList.listIterator();

		while(agentIter.hasNext()){
			anAgent = agentIter.next();
			if(anAgent.isDead){
				agentIter.remove();
				agentList.remove(anAgent);
			}
		}


	}

	/**
	 * \brief Remove an agent from the grid
	 * 
	 * Remove an agent from the grid
	 * 
	 * @param anAgent	Agent to be removed from the grid
	 */
	public void removeLocated(SpecialisedAgent anAgent) {
		if (anAgent instanceof LocatedAgent) {
			LocatedAgent aLoc = (LocatedAgent) anAgent;
			int index = getIndexedPosition(aLoc.getLocation());
			if (!Double.isNaN(index))
				_grid[index].remove(aLoc);
		}
	}

	/**
	 * \brief Update the position of the agent on the agent grid
	 * 
	 * Update the position of the agent on the agent grid. The agent may have moved due to shoving, pushing, or growth
	 * 
	 * @param anAgent	Agent that has moved and needs the location updated
	 */
	public void registerMove(LocatedAgent anAgent) {
		// Compute the theoretical index on the agentGrid
		int newIndex = getIndexedPosition(anAgent.getLocation());
		int oldIndex = anAgent.getGridIndex();

		// If gridIndex has changed, update the references
		if (isValid(anAgent.getLocation())) {
			if (newIndex != oldIndex) {
				_grid[oldIndex].remove(anAgent);
				_grid[newIndex].add(anAgent);
				anAgent.setGridIndex(newIndex);
			}
		} else {
			utils.LogFile.writeLog("Agent location is not valid -> Killed");
			//anAgent.death = "overBoard";
			anAgent.die(false);
		}
	}

	/**
	 * \brief Iterates through each square of the agent grid to calculate total biomass in that grid space, storing this on the supplied biomass grid
	 * 
	 * Iterates through each square of the agent grid to calculate total biomass in that grid space, storing this on the supplied biomass grid
	 * 
	 * @param biomassGrid	The biomass grid that will contain the total biomass in each grid square of the agent grid
	 */
	public void fitAgentMassOnGrid(SpatialGrid biomassGrid) 
	{
		for (int i = 0; i < _nTotal; i++) 
		{
			for (LocatedAgent aLoc : _grid[i].group) 
			{
				aLoc.fitMassOnGrid(biomassGrid);
			}
		}
	}

	/**
	 * \brief Iterates through each square of the agent grid to calculate total volume in that grid space, storing this on the supplied biomass grid
	 * 
	 * Iterates through each square of the agent grid to calculate total volume in that grid space, storing this on the supplied biomass grid
	 * 
	 * @param biomassGrid	The biomass grid that will contain the total volume in each grid square of the agent grid
	 */
	public void fitAgentVolumeRateOnGrid(SpatialGrid biomassGrid) {
		biomassGrid.setAllValueAt(0d);
		for (int i = 0; i < _nTotal; i++) {
			for (LocatedAgent aLoc : _grid[i].group) {
				aLoc.fitVolRateOnGrid(biomassGrid);
				;
			}
		}
	}



	/* ____________________ REPORT FILE EDITION__________________________ */

	/**
	 * \brief Summarise the agents in this container in simulation output, the form of a total biomass grid
	 * 
	 * Summarise the agents in this container in simulation output, the form of a total biomass grid
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param bufferState	The agent_state result file output buffer
	 * @param bufferSum	The agent_sum result file output buffer
	 * @throws Exception	Exception thrown if there are issues writing to these buffers
	 */
	public void writeGrids(Simulator aSim, ResultFile bufferState,
			ResultFile bufferSum) throws Exception {
		LocatedAgent aLoc;

		//sonia:chemostat
		//I've modified the refreshElement() method for the chemostat case

		// Refresh the space occupation map (0:carrier,1:biofilm or 2:liquid)
		for (int index = 0; index < _nTotal; index++) {
			_grid[index].refreshElement();
		}


		/* Build a grid of biomass concentration */

		// Set existing grid to zero
		for (int iSpecies = 0; iSpecies < aSim.speciesList.size(); iSpecies++)
			_speciesGrid[iSpecies].setAllValueAt(0);

		// Sum biomass concentrations
		for (SpecialisedAgent anA : agentList) {
			if (anA instanceof LocatedAgent) {
				aLoc = (LocatedAgent) anA;
				aLoc.fitMassOnGrid(_speciesGrid[aLoc.speciesIndex]);

			}
		}

		// now output the biomass values
		for (int iSpecies = 0; iSpecies < aSim.speciesList.size(); iSpecies++) 
		{
			_speciesGrid[iSpecies].writeReport(bufferState, bufferSum);
		}
	}

	/**
	 * \brief Write simulation output to XML files, summarising information for each species and plasmid on the agent grid
	 * 
	 * Write simulation output to XML files, summarising information for each species and plasmid on the agent grid. These files are 
	 * written out as agent_state([simulation step number]).XML and agent_sum([simulation step number]).XML
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param bufferState	The agent_state result file output buffer
	 * @param bufferSum	The agent_sum result file output buffer
	 * @throws Exception	Exception thrown if there are issues writing to these buffers
	 */
	public void writeReport(Simulator aSim, ResultFile bufferState, ResultFile bufferSum) throws Exception {

		// bvm 10.2.2009: include information about the shoving grid
		StringBuffer gridInfo = new StringBuffer();
		gridInfo.append("<grid");
		gridInfo.append(" resolution=\"").append(_res).append("\"");
		gridInfo.append(" nI=\"").append(_nI).append("\"");
		gridInfo.append(" nJ=\"").append(_nJ).append("\"");
		gridInfo.append(" nK=\"").append(_nK).append("\"");
		gridInfo.append("/>\n");
		bufferState.write(gridInfo.toString());
		bufferSum.write(gridInfo.toString());

		//sonia: 5-05-09
		//include the information about the environment status at the beginning of the agent_State and
		//agent_Sum files
		if(Simulator.isFluctEnv){
			StringBuffer envInfo = new StringBuffer();
			envInfo.append("<Environment");
			envInfo.append(" env=\"").append(FluctEnv.envStatus).append("\"");
			envInfo.append("/>\n");
			bufferSum.write(envInfo.toString());
			bufferState.write(envInfo.toString());}


		// Detail the header
		Species aSpecies;
		int spIndex, nSpecies;

		nSpecies = aSim.speciesList.size();
		StringBuffer[] speciesBuffer = new StringBuffer[nSpecies];

		// Initialise a Species markup for each present species
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			aSpecies = aSim.speciesList.get(iSpecies);
			speciesBuffer[iSpecies] = new StringBuffer();

			speciesBuffer[iSpecies].append("<species name=\"");
			speciesBuffer[iSpecies].append(aSpecies.speciesName).append("\" ");

			speciesBuffer[iSpecies].append("header=\"");
			speciesBuffer[iSpecies].append(
					aSpecies.getProgenitor().sendHeader()).append("\"");
			speciesBuffer[iSpecies].append(">\n");
		}

		// Initialise statistics (population total mass, growth-rate)
		double[] spPop = new double[nSpecies];
		double[] spMass = new double[nSpecies];
		double[] spGrowth = new double[nSpecies];


		/* <----- HGT Stats Begin ---------> */



		//sonia 26.02.2010 - Plasmid stats
		double[] spPlasmid = new double[nSpecies];
		double[] spConjugEvents = new double[nSpecies];
		int plasmidListSize=0;

		if(Simulator.multiEpi){
			plasmidListSize = aSim.plasmidList.size();
		}
		double [][][] spPlasmidTypes = new double [nSpecies][plasmidListSize][1];

		if(Simulator.multiEpi){
			for (int i = 0; i < nSpecies; i++) {
				spPlasmid[i]=0;
			}
			//sonia: 08-05-09
			for (int r=0; r< spPlasmidTypes.length; r++){
				//System.out.println("spPlasmidTypes list size is " + spPlasmidTypes.length);
				for (int c=0; c<spPlasmidTypes[r].length; c++){
					//System.out.println("spPlasmidTypes[PL INDEX] list size is " + spPlasmidTypes[r].length);
					spPlasmidTypes[r][c][0]=0;
					//System.out.println("spPlasmidTypes[r][c][0] value is " + spPlasmidTypes[r][c][0]);

				}
			}
		}



		/* <----- HGT Stats End ----> */


		for (int i = 0; i < nSpecies; i++) {
			spPop[i] = 0;
			spMass[i] = 0;
			spGrowth[i] = 0;
		}

		// Fill the agent_state file, build the state for the summary
		LocatedAgent aLoc;
		MultiEpiBac anEpiBac;

		for (SpecialisedAgent anAgent : agentList) {
			spIndex = anAgent.getSpecies().speciesIndex;
			spPop[spIndex]++;

			if (anAgent instanceof LocatedAgent & !anAgent.isDead) {
				aLoc = (LocatedAgent) anAgent;	
				spMass[spIndex] += aLoc.getTotalMass();
				spGrowth[spIndex] += aLoc.getNetGrowth();	

				/*<-------HGT Sonia Stats Begin ------> */
				//sonia:07-05-09
				if(aLoc instanceof MultiEpiBac){
					anEpiBac = (MultiEpiBac) aLoc;
					if (anEpiBac._plasmidHosted.size() != 0){
						spPlasmid[spIndex]++;
						for (int i=0; i< anEpiBac._plasmidHosted.size(); i++){
							for(int j=0; j< aSim.plasmidList.size(); j++){
								if (anEpiBac._plasmidHosted.get(i).getName().equals(aSim.plasmidList.get(j))){
									int	plIndex = aSim.plasmidList.indexOf(aSim.plasmidList.get(j));
									spPlasmidTypes[spIndex][plIndex][0]++;		
									//System.out.println("plasmid count is  " + spPlasmidTypes[spIndex][plIndex][0]);
								}
							}
						}
					}

					if(anEpiBac.plasmidVector.size()!= 0){
						spConjugEvents[spIndex]= anEpiBac.plasmidVector.size() + spConjugEvents[spIndex];
					}

				}else{
					spPlasmid[spIndex]=0;

				}
				/* <------- HGT Sonia Stats End ------> */

				speciesBuffer[spIndex].append(aLoc.writeOutput()+";\n");

			}
		}
		
		
		StringBuffer text;
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			text = new StringBuffer();
			aSpecies = aSim.speciesList.get(iSpecies);
			text.append("<species name=\"");
			text.append(aSpecies.speciesName).append("\" ");
			text.append("header=\"");
			bufferSum.write(text.toString());

			/*<----HGT Sonia begin---->*/
			if(Simulator.multiEpi){
				text= new StringBuffer("population,mass,growthRate,plasmidBearing");

				if(aSim.plasmidList.size()>0){
					for (int iPl = 0; iPl< aSim.plasmidList.size(); iPl++){
						String plName = aSim.plasmidList.get(iPl);
						text.append(",");
						text.append(plName);
					}
				}			
				text.append(",conjCount");
				/*<----HGT Sonia end---->*/

			}else{
				text= new StringBuffer("population,mass,growthRate");
			}

			text.append("\" ");
			text.append(">\n");

			bufferSum.write(text.toString());

			//writing stats on population, mass and growth
			text = new StringBuffer ("");
			text.append(spPop[iSpecies] + "," + spMass[iSpecies] + ","
					+ spGrowth[iSpecies] );

			/*<----HGT Sonia begin---->*/
			if(Simulator.multiEpi){
				text.append("," + spPlasmid[iSpecies]);
				for (int c=0; c<spPlasmidTypes[iSpecies].length; c++){
					double plCount = spPlasmidTypes[iSpecies][c][0];
					text.append("," + plCount);			
				}
				text.append("," + spConjugEvents[iSpecies] + ";");
			}
			/*<----HGT Sonia end---->*/

			text.append("</species>\n");
			bufferSum.write(text.toString());
		}

		//brian
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			speciesBuffer[iSpecies].append("</species>\n");
			bufferState.write(speciesBuffer[iSpecies].toString());
		}	

	}


	/**
	 * \brief Output information on agents that have been removed from the simulation. Will enable investigators to find out why this is the case
	 * 
	 * Output information on agents that have been removed from the simulation. Will enable investigators to find out why this is the case
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param bufferStateDeath	The agent_state_death result file output buffer
	 * @param bufferSumDeath	The agent_sum_death result file output buffer
	 * @throws Exception	Exception thrown if there are issues writing to these buffers
	 */
	public void writeReportDeath(Simulator aSim, ResultFile bufferStateDeath, ResultFile bufferSumDeath) throws Exception 
	{
		
		// STATUS MESSAGE TO THE SCREEN, IF NOT IN QUIET MODE
		if (!Idynomics.quietMode)
			System.out.println("size of agentToKill list at beginning of the writeReportDeath:  " + _agentToKill.size());

		// DEATH FILE FORMAT (FOR BOTH STATE AND SUM)
		// FIRST SET OF TAGS LISTS THE SIMULATION ITERATION, AND THE TIME THIS REPRESENTS
		// THIS IS WRITTEN BY THE METHOD THAT CALLS THIS, NOT WRITTEN HERE.
		// THUS THIS STARTS BY WRITING INFORMATION ABOUT THE SHOVING GRID
				
		StringBuffer gridInfo = new StringBuffer();
		// Generate the grid info:
		gridInfo = writeGridInformation(gridInfo);
		// Write this to both the agent death and removed biomass files
		bufferSumDeath.write(gridInfo.toString());
		bufferStateDeath.write(gridInfo.toString());

		// IF WE ARE DEALING WITH A FLUCTENV SIMULATION, THERE ARE ADDITIONAL ENVIRONMENTAL DATA TO WRITE TO THE FILES
		if(Simulator.isFluctEnv)
		{
			StringBuffer envInfo = writeFluctEnvInformation();
			bufferSumDeath.write(envInfo.toString());
			bufferStateDeath.write(envInfo.toString());
		}

		// NOW WE'RE GOING TO CREATE A HEADER FOR EACH OF THE SPECIES IN THE SIMULATION
		// HEADER WILL CONTAIN SPECIES NAME, FOLLOWED BY THE NAMES OF EACH RESULT COLUMN
		Species aSpecies;
		int spIndex, nSpecies;

		// Get the number of species in this simulation
		nSpecies = aSim.speciesList.size();
		// Set up a buffer to hold information for each agent of these species
		StringBuffer[] speciesBuffer = new StringBuffer[nSpecies];

		// Now to iterate through each species and create a markup for each
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) 
		{
			aSpecies = aSim.speciesList.get(iSpecies);
			speciesBuffer[iSpecies] = new StringBuffer();

			// First part, the species name
			speciesBuffer[iSpecies].append("<species name=\"");
			speciesBuffer[iSpecies].append(aSpecies.speciesName).append("\" ");

			// Now create the header. Note that this information comes from the species itself - its part of the agent and the inherited
			// classes. So if we are going to change that, it gets changed there
			speciesBuffer[iSpecies].append("header=\"");
			speciesBuffer[iSpecies].append(
					aSpecies.getProgenitor().sendHeader());
			
			// Now we can add additional fields if these are not declared in the agent itself. In this case, we're going to add a reason
			// for that agents death
			speciesBuffer[iSpecies].append(",death");
			speciesBuffer[iSpecies].append("\"");
			speciesBuffer[iSpecies].append(">\n");
		}

		// REMEMBER WE ARE CREATING TWO STATS FILES AT THE SAME TIME HERE - ONE CONTAINING ALL THE AGENTS THAT ARE DEAD, AND ONE CONTAINING
		// STATISTICS FOR DEATH IN THIS OUTPUT ITERATION. HERE WE INITIALISE STORAGE FOR STATISTICS (POPULATION, TOTAL MASS, GROWTH RATE)
		// TO BE CONTAINED IN THE AGENTSUM_DEATH FILE
		double[] spPop = new double[nSpecies];
		double[] spMass = new double[nSpecies];
		double[] spGrowth = new double[nSpecies];
		// INITIALISE THESE ARRAYS
		for (int i = 0; i < nSpecies; i++) 
		{
			spPop[i] = 0;
			spMass[i] = 0;
			spGrowth[i] = 0;
		}

		// FOR PLASMID SIMULATIONS< THERE ARE EXTRA STATISTICS. THIS SECTION DEALS WITH THIS
		/* <----- HGT Stats Begin ---------> */
		double[] spPlasmid = new double[nSpecies];
		double[] spConjugEvents = new double[nSpecies];
		int plasmidListSize=0;

		if(Simulator.multiEpi)
		{
			plasmidListSize = aSim.plasmidList.size();
		}
		double [][][] spPlasmidTypes = new double [nSpecies][plasmidListSize][1];

		if(Simulator.multiEpi)
		{
			for (int i = 0; i < nSpecies; i++) 
			{
				spPlasmid[i]=0;
			}
			//sonia: 08-05-09
			for (int r=0; r< spPlasmidTypes.length; r++)
			{
				for (int c=0; c<spPlasmidTypes[r].length; c++)
				{
					spPlasmidTypes[r][c][0]=0;
				}
			}
		}	
		/* <----- HGT Stats End ----> */	

		// NOW TO GET BACK TO THE NORMAL SIMULATION CASE, LETS BUILD THE AGENT STATE DEATH FILE AND STATISTICS
		LocatedAgent aLoc;
		MultiEpiBac anEpiBac;
		
		// FOR EACH AGENT IN THE _AGENTTOKILL LIST
		for (SpecialisedAgent anAgent : _agentToKill) 
		{
			spIndex = anAgent.getSpecies().speciesIndex;
			spPop[spIndex]++;

			if (anAgent instanceof LocatedAgent) 
			{
				aLoc = (LocatedAgent) anAgent;	
				// Increment the total mass and total growth rate for this agent
				spMass[spIndex] += aLoc.getTotalMass();
				spGrowth[spIndex] += aLoc.getNetGrowth();	

				// Now examine the statistics for a plasmid simulation and increase these where necessary
				/* <----- HGT Sonia Stats Start ----> */
				if(aLoc instanceof MultiEpiBac)
				{
					anEpiBac = (MultiEpiBac) aLoc;
					if (anEpiBac._plasmidHosted.size() != 0)
					{
						spPlasmid[spIndex]++;
						for (int i=0; i< anEpiBac._plasmidHosted.size(); i++)
						{
							for(int j=0; j< aSim.plasmidList.size(); j++)
							{
								if (anEpiBac._plasmidHosted.get(i).getName().equals(aSim.plasmidList.get(j)))
								{
									int	plIndex = aSim.plasmidList.indexOf(aSim.plasmidList.get(j));
									spPlasmidTypes[spIndex][plIndex][0]++;		
								}
							}
						}
					}

					if(anEpiBac.plasmidVector.size()!= 0){
						spConjugEvents[spIndex]= anEpiBac.plasmidVector.size() + spConjugEvents[spIndex];
					}

				}
				else
				{
					spPlasmid[spIndex]=0;
				}
				/* <----- HGT Sonia Stats End ----> */

				// NOW THE STATS ARE DEALT WITH, WRITE THE INFORMATION FOR THIS AGENT OF THIS SPECIES TO THE BUFFER
				// NOTE THAT THE CODE TO DO THIS IS CONTAINED IN THE AGENT ITSELF
				speciesBuffer[spIndex].append(aLoc.writeOutput());
				
				// NOW ADD ANY ADDITIONAL INFORMATION THAT WAS ADDED IN THE HEADER PREVIOUSLY 
				// IN THIS CASE REASON FOR DEATH
				speciesBuffer[spIndex].append("," + aLoc.death + ";\n");

			}
		}
		
		

		// NOW WE HAVE THE STATS AND HAVE GONE THROUGH EACH AGENT BEING KILLED, WE CAN CREATE THE AGENT_SUM_DEATH FILE
		StringBuffer text;
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) 
		{
			text = new StringBuffer();
			aSpecies = aSim.speciesList.get(iSpecies);
			text.append("<species name=\"");
			text.append(aSpecies.speciesName).append("\" ");
			text.append("header=\"");
			bufferSumDeath.write(text.toString());

			text= new StringBuffer("population,mass,growthRate");

			//<--- Multi HGT starts---->
			if (Simulator.multiEpi){			
				if(aSim.plasmidList.size()>0){
					text.append(",plasmidBearing");
				}

				if(aSim.plasmidList.size()>0){
					for (int iPl = 0; iPl< aSim.plasmidList.size(); iPl++){
						String plName = aSim.plasmidList.get(iPl);
						text.append(",");
						text.append(plName);
					}
				}

				if(aSim.plasmidList.size()>0){
					text.append(",conjCount");
				}
			}
			//<--- Multi HGT ends---->

			text.append("\" ");
			text.append(">\n");

			bufferSumDeath.write(text.toString());

			text = new StringBuffer ("");
			text.append(spPop[iSpecies] + "," + spMass[iSpecies] + ","
					+ spGrowth[iSpecies]);

			//<--- Multi HGT starts---->
			if (Simulator.multiEpi){
				if(aSim.plasmidList.size()>0){
					text.append("," + spPlasmid[iSpecies]);
				}

				if(aSim.plasmidList.size()>0){
					for (int c=0; c<spPlasmidTypes[iSpecies].length; c++){
						double plCount = spPlasmidTypes[iSpecies][c][0];
						text.append("," + plCount);			
					}
				}

				if(aSim.plasmidList.size()>0){
					text.append("," + spConjugEvents[iSpecies]);
					text.append("</species>\n");
				}
				//<--- Multi HGT ends---->
			}else{
				//text.append(";");
				text.append("</species>\n");
			}

			bufferSumDeath.write(text.toString());
		}
		
		// NOW CAN ITERATE THROUGH THE SPECIES AND WRITE THE AGENT_STATE_DEATH FILE FOR EACH AGENT OF THAT SPECIES
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) 
		{
			speciesBuffer[iSpecies].append("</species>\n");
			bufferStateDeath.write(speciesBuffer[iSpecies].toString());
		}	

		


	}

	/**
	 * \brief Writes information on the shoving grid to the output string specified
	 * 
	 * Writes information on the shoving grid to the output string specified
	 * 
	 * @param outputString	String that will be output in the result file
	 * @return	String with information about the grid (resolution and size in I,J,K directions)
	 */
	public StringBuffer writeGridInformation(StringBuffer outputString)
	{
		outputString.append("<grid");
		outputString.append(" resolution=\"").append(_res).append("\"");
		outputString.append(" nI=\"").append(_nI).append("\"");
		outputString.append(" nJ=\"").append(_nJ).append("\"");
		outputString.append(" nK=\"").append(_nK).append("\"");
		outputString.append("/>\n");
		
		return outputString;
	}
	
	/**
	 * \brief Writes information on the environment to the result files if this is a fluctenv simulation
	 * 
	 * Writes information on the environment to the result files if this is a fluctenv simulation
	 * 
	 * @return	Stringbuffer with information about the environment for writing to the output file
	 */
	public StringBuffer writeFluctEnvInformation()
	{
		StringBuffer envInfo = new StringBuffer();
		envInfo.append("<Environment");
		envInfo.append(" env=\"").append(FluctEnv.envStatus).append("\"");
		envInfo.append("/>\n");
		
		return envInfo;
	}


	/*public void preprintLevelSet() {
		// Build the matrix of erosion time
		// _levelset.refreshFromBorder();

		for (int index = 0; index < _nTotal; index++) {
			_grid[index].printLevelSet(_erosionGrid);
		}
	}*/


	/*public void printLevelSet(ResultFile bufferState) throws Exception {

		// Edit the markup for the solute grid
		StringBuffer value = new StringBuffer();
		value.append("<solute name=\"").append("erosionSpeed");
		value.append("\" unit=\"hr");
		value.append("\" resolution=\"").append(_res);
		value.append("\" nI=\"").append(_nI);
		value.append("\" nJ=\"").append(_nJ);
		value.append("\" nK=\"").append(_nK);
		value.append("\">\n");

		// Write the markup in the file
		bufferState.write(value.toString());

		// Fill the mark-up
		if (_nK == 1) {
			// We have a 2D grid
			for (int i = 0; i < _nI+2; i++) {
				for (int j = 0; j < _nJ+2; j++) {
					bufferState.write(Arrays.toString(_erosionGrid[i][j]));
					bufferState.write(";\n");
				}
			}
		} else {
			for (int i = 0; i < _nI+2; i++) {
				for (int j = 0; j < _nJ+2; j++) {
					bufferState.write(Arrays.toString(_erosionGrid[i][j]));
					bufferState.write(";\n");
				}
			}
		}

		// Close the mark-up
		bufferState.write("\n</solute>\n");
	}*/

	/**
	 * \brief Checks the grid size, dependent on whether this is a biofilm or chemostat simulation, and sets suitable dimensions
	 * 
	 * Checks the grid size, dependent on whether this is a biofilm or chemostat simulation, and sets suitable dimensions
	 * 
	 * @param aSimulator	The simulation object used to simulate the conditions specified in the protocol file
	 * @param root	The agentGrid markup from the XML file
	 */
	public void checkGridSize(Simulator aSimulator, XMLParser root) 
	{
		// Read in the grid resolution from the XML protocol file
		_res = root.getParamDbl("resolution");

		if(Simulator.isChemostat)
		{
			//set the resolution to the resolution of the domain
			_res = domain._resolution;
			//do not correct the grid size
			//set the number of grid elements to 1
			_nTotal = 1;
		}
		else
		{
			// Eventually correct grid size
			
			_nI = (int) Math.ceil(domain.length_X / _res);
			_res = domain.length_X / _nI;
			_nJ = (int) Math.ceil(domain.length_Y / _res);
			
			// Now determine if dealing with a 3D environment
			if (domain.is3D()) 
			{
				_nK = (int) Math.ceil(domain.length_Z / _res);
				is3D = true;
				
				// KA 050613 - changed this here to create the _nTotal dependent on whether we are in 3D or not - little need for
				// Z padding if we're not in 3D
				// KA 270713 - bug found that creates different results between 1.1 and 1.2 - although this has been found and fixed
				// this appears padding related, and thus 1.2 has returned to have padding for 2D simulations. This will be addressed 
				// properly at a later version
				
				// Calculate the number of grid elements
				_nTotal = (_nI + 2) * (_nJ + 2) * (_nK + 2);
				
			} 
			else 
			{
				_nK = 1;
				is3D = false;
				
				// KA 270713 - bug found that creates different results between 1.1 and 1.2 - although this has been found and fixed
				// this appears padding related, and thus 1.2 has returned to have padding for 2D simulations. This will be addressed 
				// properly at a later version

				// KA 260713 - to check for differences between 1.1 and 1.2, have put padding back into 2D sim
				_nTotal = (_nI + 2) * (_nJ + 2) * (_nK + 2);
				
				// Code for grid without padding - turn on again when this is investigated further
				// Now work out _Total without the addition of padding in the Z, which is not necessary and should save space
				//_nTotal  = (_nI + 2) * (_nJ + 2);
				
				
			}
		}

		
	}

	/**
	 * \brief Creates the shoving grid, padded of LocatedGroups
	 * 
	 * Create a vectorized array of spatial groups, build their neighbourhood and store it 
	 * Note :Shoving grid is a padded of LocatedGroups
	 * 
	 * @param aSimulator	The current simulation object that is recreating the conditions specified in the protocol file
	 */
	public void createShovGrid(Simulator aSimulator) 
	{
		_grid = new LocatedGroup[_nTotal];	
		for (int index = 0; index < _nTotal; index++)
			_grid[index] = new LocatedGroup(index, this, aSimulator);

		for (int index = 0; index < _nTotal; index++)
			_grid[index].init();


	}

	/**
	 * \brief Creates the output species and erosion grids
	 * 
	 * Creates the output species and erosion grids. There is one erosion grid, but one grid for each species in the simulation
	 * 
	 * @param aSim	The current simulation object that is recreating the conditions specified in the protocol file
	 */
	public void createOutputGrid(Simulator aSim) 
	{
		_erosionGrid = new double[_nI + 2][_nJ + 2][_nK + 2];
		_speciesGrid = new SpatialGrid[aSim.speciesList.size()];
		
		// Create a grid for each species in the simulation
		for (int iSpecies = 0; iSpecies < aSim.speciesDic.size(); iSpecies++)
			_speciesGrid[iSpecies]=domain.createGrid(aSim.speciesDic.get(iSpecies), 0);
	}

	/* _____________________ EROSION & DETACHMENT ___________________________ */

	/**
	 * \brief Perform connected volume filtration (connected to bottom) to determine which agents should be marked for sloughing
	 * 
	 *  Perform connected volume filtration (connected to bottom) to determine which agents should be marked for sloughing
	 */
	protected void markForSloughing() {
		// 
		// (cvf is true for connected elements, and false for non-connected)
		boolean[] cvf = (new ConnectedVolume(_nI, _nJ, _nK)).computeCvf(_grid);

		int nRemoved = 0;
		double mRemoved = 0;

		// mark as detachable all particles in non-valid map positions
		for (int index = 0; index < _nTotal; index++) {
			// if (_grid[index].status==1&&!cvf[index]) {
			if (!cvf[index]) {
				// if it's not connected, remove the agents
				if (_grid[index].totalMass > 0) {
					nRemoved += _grid[index].group.size();
					mRemoved += _grid[index].totalMass;

					_grid[index].killAll();
				}
			}
		}

		LogFile.writeLog("Sloughing " + nRemoved + " ("
				+ ExtraMath.toString(mRemoved, false) + " fg)");
	}

	/**
	 * \brief Find the border points which go through a process of erosion (erosion and sloughing have changed the configuration)
	 * 
	 * Find the border points which go through a process of erosion (erosion and sloughing have changed the configuration)
	 * 
	 */
	public void shrinkOnBorder() 
	{
		_levelset.refreshBorder(true, mySim);

		double mass = 0;
		double tallyVariable = 0;
		int nDetach = 0;
		int index;
		double ratio;

		for (LocatedGroup aBorderElement : _levelset.getBorder()) 
		{
			index = aBorderElement.gridIndex;

			ratio = SimTimer.getCurrentTimeStep() / _grid[index].erosionTime;
			ratio = Math.max(0, Math.min(ratio, 1));
			tallyVariable = aBorderElement.totalMass * ratio;


			for (LocatedAgent aLoc : _grid[index].group) {
				mass += aLoc.getTotalMass() * ratio;
				for (int iComp = 0; iComp < aLoc.particleMass.length; iComp++)
					aLoc.particleMass[iComp] *= 1 - ratio;

				aLoc.updateSize();
				if (aLoc.willDie()) {
					mass += aLoc.getTotalMass();
					aLoc.die(false);
					aLoc.death = "detachment";

					nDetach++;
				}
			}
			if (mass > tallyVariable)
				continue;
		}

		LogFile.writeLog("Eroding " + nDetach + " ("
				+ ExtraMath.toString(mass, true) + "/"
				+ ExtraMath.toString(tallyVariable, true) + " fg)");
	}
	

	/**
	 * \brief Calculate the detachment priority
	 * 
	 * @param i
	 * @param location
	 * @return
	 */
	private double detFunction(double i, double location) {
		if (i<1){ return ExtraMath.sq(_res - (location % _res));}
		else {return ExtraMath.sq(location % _res);}
	}

	/**
	 * \brief Performs the discrete removal of agents over the whole _close group together, instead of considering each Located Group seperately. 
	 * 
	 * Performs the discrete removal of agents over the whole _close group together, instead of considering each Located Group 
	 * seperately. Note also that here tallyVariable is cumulative: whatever remains rolls over to the next time-step. This eliminates 
	 * any timestep issues.
	 * 
	 * @param agentGrid	The agent grid
	 */
	public void removeOnBorder(AgentContainer agentGrid) 
	{
		// Find the border points (erosion and sloughing have changed the
		// configuration)
		_levelset.refreshBorder(true, mySim);

		double mass = 0; // counter of mass so far removed
		int nDetach = 0; // counter of agents removed
		LinkedList<LocatedAgent> detGroup = new LinkedList<LocatedAgent>(); // List of agents to consider for removal
		LocatedAgent detLoc; // the member of detGroup with the highest detPriority

		// For groups on _close list
		for (LocatedGroup aBorderElement : _levelset.getBorder()) {

			// tally up tallyVariable, the approximate amount of mass to remove, by the ratio variable of each border element
			aBorderElement.ratio = SimTimer.getCurrentTimeStep() / _grid[aBorderElement.gridIndex].erosionTime;
			// i.e. ratio = (currentTimeStep * detachmentSpeed * numberFreeNeighbours)/(resolution)
			aBorderElement.ratio = Math.min(aBorderElement.ratio, 1);
			tallyVariable += aBorderElement.totalMass*aBorderElement.ratio;

			// and add them to detGroup
			for (LocatedAgent aLoc:aBorderElement.group) detGroup.add(aLoc);

		} // end of: for (LocatedGroup aBorderElement : _levelset.getBorder())

		if (tallyVariable>Collections.min(detGroup, new LocatedAgent.totalMassComparator()).getTotalMass())  {

			// For groups on _close list
			for (LocatedGroup aBorderElement : _levelset.getBorder()) {
				// calculate detPriority for all agents in the _close list
				calcDetPriority(agentGrid, aBorderElement, aBorderElement.ratio);
			}
			
			detLoc = Collections.max(detGroup, new LocatedAgent.detPriorityComparator());
			while (detLoc.getTotalMass()<tallyVariable && detGroup.size()>0) {
				// Remove agents one by one, in order of decreasing detPriority

				detLoc = Collections.max(detGroup, new LocatedAgent.detPriorityComparator());
				mass += detLoc.getTotalMass();
				tallyVariable -= detLoc.getTotalMass();
				detLoc.die(false);
				detLoc.death = "detachment";
				nDetach++;
				detGroup.remove(detLoc);


			} // end of if (2*tallyVariable>Collections.min(aBorderElement.group, new LocatedAgent....)
		}

		System.out.println("******************************REMOVE ON BORDER**************************************");
		LogFile.writeLog("Eroding " + nDetach + " ("
				+ ExtraMath.toString(mass, true) + "/"
				+ ExtraMath.toString(tallyVariable, true) + " fg) from "
				+ _levelset.getBorder().size() +" elements.");
	}


	/**
	 * \brief Calculate detachment priorties
	 * 
	 * Calculate detachment priorties
	 * 
	 * @param agentGrid	This agent grid of located agents
	 * @param aBorderElement	Group of located agents that are on the border
	 * @param ratio	Ratio value set for that LocatedGroup
	 */
	public void calcDetPriority(AgentContainer agentGrid, LocatedGroup aBorderElement, double ratio)
	{
		int i = 0;

		// Reset all detPriority values to zero
		for (LocatedAgent aLoc:aBorderElement.group) aLoc.detPriority=0;

		// For each free neighbour run through the agents in our border element, 
		// adding the square of the agent's proximity to that neighbour
		// (i.e. distance from opposite neighbour) to its detPriority variable
		for (i=0;i<3;i+=2){
			// x-side
			if (aBorderElement.nbhGroup[i][1][1].status==2) {
				// LogFile.writeLog(aBorderElement.nbhGroup[i][1][1].dc+"is free");
				for (LocatedAgent aLoc:aBorderElement.group) {
					aLoc.detPriority += detFunction(i,aLoc.getLocation().x);
					// LogFile.writeLog("Agent: "+aLoc.sendName()+", at "+aLoc.getLocation().x+", detPriority: "+aLoc.detPriority);
				}
			}
			// y-side & z-side (double y-side if 2D)
			if (agentGrid.is3D) {
				if (aBorderElement.nbhGroup[1][i][1].status==2) {
					// LogFile.writeLog(aBorderElement.nbhGroup[1][i][1].dc+"is free");
					for (LocatedAgent aLoc:aBorderElement.group) {
						aLoc.detPriority += detFunction(i,aLoc.getLocation().y);
						// LogFile.writeLog("Agent: "+aLoc.sendName()+", at "+aLoc.getLocation().y+", detPriority: "+aLoc.detPriority);
					}
				}
				if (aBorderElement.nbhGroup[1][1][i].status==2) {
					// LogFile.writeLog(aBorderElement.nbhGroup[1][1][i].dc+"is free");
					for (LocatedAgent aLoc:aBorderElement.group) {
						aLoc.detPriority += detFunction(i,aLoc.getLocation().z);
						// LogFile.writeLog("Agent: "+aLoc.sendName()+", at "+aLoc.getLocation().z+", detPriority: "+aLoc.detPriority);
					}
				}
			} else {
				if (aBorderElement.nbhGroup[1][i][1].status==2) {
					// LogFile.writeLog(aBorderElement.nbhGroup[1][i][1].dc+"is free");
					for (LocatedAgent aLoc:aBorderElement.group) {
						aLoc.detPriority += 2*detFunction(i,aLoc.getLocation().y);
						// LogFile.writeLog("Agent: "+aLoc.sendName()+", at "+aLoc.getLocation().y+", detPriority: "+aLoc.detPriority);
					} // end of: for (LocatedAgent aLoc:aBorderElement.group) {
				} // end of: if (aBorderElement.nbhGroup[1][i][1].status==2) {
			} // end of: if (agentGrid.is3D) {...}else{
		} // end of: for (i=0;i<3;i+=2){

		// weight the detPriority of each agent by its Located Group ratio
		for (LocatedAgent aLoc:aBorderElement.group) aLoc.detPriority *= ratio;
	}

	/* __________________________ GET & SET _________________________________ */

	/**
	 * \brief Check that the coordinates expressed in a DiscreteVector are defined on this grid (the grid is padded but dc uses shifted coordinates)
	 * 
	 * Check that these coordinates expressed in a DiscreteVector are defined on this grid (the grid is padded but dc uses shifted coordinates)
	 * 
	 * @return Boolean stating whether or not the DiscreteVector is valid
	 */
	public boolean isValid(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i >= 0) & (dC.i < _nI);
		out &= (dC.j >= 0) & (dC.j < _nJ);
		out &= (dC.k >= 0) & (dC.k < _nK);
		return out;
	}

	/**
	 * \brief Check that the coordinates expressed in a ContinuousVector are defined on this grid (the grid is padded but dc uses shifted coordinates)
	 * 
	 * Check that these coordinates expressed in a ContinuousVector are defined on this grid (the grid is padded but dc uses shifted coordinates)
	 * 
	 * @return Boolean stating whether or not the ContinuousVector are valid
	 */
	public boolean isValid(ContinuousVector cC) {
		return isValid(getGridPosition(cC));
	}



	/**
	 * \brief Find the voxel a continuous coordinate lies in and return the index
	 * 
	 * Find the voxel a continuous coordinate lies in and return the index
	 * 
	 * @param cc	the continuous coordinate to find the index for
	 * @return index on the 1D (vectorized) array
	 */
	public int getIndexedPosition(ContinuousVector cc) 
	{

		//sonia:chemostat
		// to guarantee that the agents are effectively removed
		if(Simulator.isChemostat){
			return 0;

		}
		else
		{
			// KA 050613 - changed when examining why Z was created with padding when using 2D - now calculates dependent on dimension
			int i = (int) Math.floor(cc.x / _res) + 1;
			int j = (int) Math.floor(cc.y / _res) + 1;
			//System.out.println("I: "+i+" J: "+j+" _nI "+_nI+" NJ"+_nJ);
			
			// KA 260713 - to check for differences between 1.1 and 1.2, have put padding back into 2D sim
			int k = (int) Math.floor(cc.z / _res) + 1;
			return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
			
			// WITH BUG BETWEEN 1.1 AND 1.2 SUGGESTING THIS NEEDS FURTHER EXAMINING, PADDING PUT BACK IN TO 2D SIMULATIONS
			// THIS CODE CAN BE RESTORED LATER WHEN THIS IS REEXAMINED
			//if(is3D)
			//{	
				//int k = (int) Math.floor(cc.z / _res) + 1;
				//return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
			//}
			//else
			//{
				//return i + j * (_nI + 2) + (_nI + 2);
			//}
		}
	}

	/**
	 * \brief Find the voxel a discrete vector coordinate lies in and return the index
	 * 
	 * Find the voxel a discrete vector coordinate lies in and return the index
	 * 
	 * @param dc	the discrete vector coordinate to find the index for
	 * @return index on the 1D (vectorized) array
	 */
	public int getIndexedPosition(DiscreteVector dc) {
		int i = dc.i + 1;
		int j = dc.j + 1;
		int k = dc.k + 1;

		return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
	}

	/**
	 * \brief Takes a ContinuousVector as input (in microns) and returns a DiscreteVector expressing its location on the grid (in voxels) 
	 * 
	 * Takes a ContinuousVector as input (in microns) and returns a DiscreteVector expressing its location on the grid (in voxels)
	 * 
	 * @param cC	ContinuousVector that expresses an agent location (in microns)
	 * @return	A discrete vector expressing the agents position in terms of the agent grid (voxels)
	 */
	public DiscreteVector getGridPosition(ContinuousVector cC) {
		int i = (int) Math.floor(cC.x / _res);
		int j = (int) Math.floor(cC.y / _res);
		int k = (int) Math.floor(cC.z / _res);

		return new DiscreteVector(i, j, k);
	}

	/**
	 * \brief Takes an voxel integer index and returns a DiscreteVector containing the X,Y, and Z coordinates of that voxel
	 * 
	 * Takes an voxel integer index and returns a DiscreteVector containing the X,Y, and Z coordinates of that voxel
	 *  
	 * @param index	Integer index specifying a voxel grid space on the agent grid
	 * @return	A discrete vector of the coordinates of this grid
	 */
	public DiscreteVector getGridPosition(int index) {
		int k = (int) Math.floor(index / (_nI + 2) / (_nJ + 2));
		int j = (int) Math.floor((index - k * ((_nI + 2) * (_nJ + 2)))
				/ (_nI + 2));
		int i = (int) Math.floor((index - (k * (_nI + 2) * (_nJ + 2)) - j
				* (_nI + 2)));

		return new DiscreteVector(i - 1, j - 1, k - 1);
	}

	/**
	 * \brief Takes an voxel integer index and returns a ContinuousVector containing the X,Y, and Z coordinates of the centre of that voxel
	 * 
	 * Takes an voxel integer index and returns a ContinuousVector containing the X,Y, and Z coordinates of the centre of that voxel
	 *  
	 * @param index	Integer index specifying a voxel grid space on the agent grid
	 * @return	A continuous vector of the coordinates of this grid
	 */
	public ContinuousVector getGridLocation(int index) {
		int k = (int) Math.floor(index / (_nI + 2) / (_nJ + 2));
		int j = (int) Math.floor((index - k * ((_nI + 2) * (_nJ + 2)))
				/ (_nI + 2));
		int i = (int) Math.floor((index - (k * (_nI + 2) * (_nJ + 2)) - j
				* (_nI + 2)));

		return new ContinuousVector((i - .5) * _res, (j - .5) * _res,
				(k - .5) * _res);
	}



	/**
	 * \brief Return the agent time step
	 * 
	 * Return the agent time step
	 * 
	 * @return	Double value stating the agent time step
	 */
	public double getAgentTimeStep() {
		return AGENTTIMESTEP;
	}

	/**
	 * \brief Return the resolution of the grid
	 * 
	 * Return the resolution of the grid
	 * 
	 * @return	Resolution of the agent grid
	 */
	public double getResolution() {
		return _res;
	}

	/**
	 * \brief Return the dimensions of the grid (nI,nJ,nK) as an array
	 * 
	 * Return the dimensions of the grid (nI,nJ,nK) as an array
	 * 
	 * @return	An array containing the dimensions of the grid
	 */
	public int[] getGridDescription() {
		int[] out = { _nI, _nJ, _nK };
		return out;
	}
	
	/**
	 * \brief Return the number of grid cells in the J direction
	 * 
	 * Return the number of grid cells in the J direction
	 * 
	 * @return	Number of grid cells in the J direction
	 */
	public int get_nJ()
	{
		return _nJ;
	}
	
	/**
	 * \brief Return the number of grid cells in the K direction
	 * 
	 * Return the number of grid cells in the K direction
	 * 
	 * @return	Number of grid cells in the K direction
	 */
	public int get_nK()
	{
		return _nK;
	}

	/**
	 * \brief Return the shoving grid
	 * 
	 * Return the shoving grid
	 * 
	 * @return	LocatedGroup containing the calculated shoving information
	 */
	
	public LocatedGroup[] getShovingGrid() {
		return _grid;
	}

	/**
	 * \brief Return the levelset used for modelling detachment
	 * 
	 * Return the levelset used for modelling detachment
	 * 
	 * @return	Levelset used for modelling detachment
	 */
	public LevelSet getLevelSet() {
		return _levelset;
	}
	
	/**
	 * \brief Return the status of a given grid voxel, to determine if the voxel is within a biofilm or not
	 * 
	 * Return the status of a given grid voxel, to determine if the voxel is within a biofilm or not
	 * 
	 * @param gridVoxel	Integer of the grid voxel whos status is being queried
	 * @return	Integer showing the status of that grid voxel
	 */
	public int getVoxelStatus(int gridVoxel)
	{
		return _grid[gridVoxel].status;
	}
	
	/**
	 * \brief Return the located group of agents in an agent grid voxel
	 * 
	 * @param gridVoxel	Integer of the grid voxel for which the located group is to be returned
	 * @return	LocatedGroup containing the agents within that grid voxel
	 */
	public LocatedGroup returnGroupInVoxel(int gridVoxel)
	{
		return _grid[gridVoxel];
	}

}
