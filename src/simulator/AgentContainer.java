/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * 
 *  Stores all the agents, call them and manages shoving/erosion of located agents
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Doetsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sonia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
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
import simulator.SpatialGrid;

import utils.ResultFile;
import utils.XMLParser;
import utils.LogFile;
import utils.ExtraMath;

public class AgentContainer {

	/* __________________ Properties __________________________ */

	public Domain domain;
	public Simulator mySim;

	// Container for all agents (even the non located ones)

	public LinkedList<SpecialisedAgent> agentList;
	public ListIterator<SpecialisedAgent> agentIter;

	// Temporary containers used to store agents who will be added or removed
	//sonia 27.04.2010
	//changed visibility to public so that it can be accessed from LocatedGroup in killAll()
	public LinkedList<SpecialisedAgent> _agentToKill = new LinkedList<SpecialisedAgent>();


	private SpatialGrid[] _speciesGrid;

	// Definition of the spatial grid ____________________________ */
	private int _nI, _nJ, _nK, _nTotal;
	private double _res;
	private LocatedGroup[] _grid;
	protected double[][][] _erosionGrid;
	public boolean is3D;

	// Parameters of the shoving algorithm _________________________________ */
	private final double SHOV_FRACTION;

	public final double AGENTTIMESTEP;
	private final int MAXITER;
	private final boolean MUTUAL;
	// boolean to select erosion method: true is shrinkOnBorder, false is removeOnBorder
	private final boolean EROSIONMETHOD;
	private final boolean DOSLOUGHING;
	private final int maxPopLimit;
	private SoluteGrid _pressure;
	
	// tally of mass to be removed in removeOnBorder, or cells in chemostat dilution (see agentsFlushedAway)
	double tallyVariable = 0; 

	// Current number and maximal number of shoving iterations
	int shovIter, shovLimit, maxShoveIter, nMoved;
	double tMoved, deltaMove;

	private LevelSet _levelset;

	//sonia 23.11.09
	//private double agentsToDilute;
	private double Dfactor;

	/* _______________________ CONSTRUCTOR __________________________________ */

	/**
	 * XML protocol file-based constructor
	 */
	public AgentContainer(Simulator aSimulator, XMLParser root, double agentTimeStep) throws Exception {
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

		//sonia 20/01/2011 - now using getParamTime() function that correctly converts time units
		//double value = root.getParamDbl("agentTimeStep");

		double value = agentTimeStep;


		if (Double.isNaN(value)){
			AGENTTIMESTEP = SimTimer.getCurrentTimeStep();
			LogFile.writeLog("Using global timestep of "+AGENTTIMESTEP+" for agentTimeStep");
		}else{
			AGENTTIMESTEP = value;
			if (AGENTTIMESTEP > SimTimer.getCurrentTimeStep()) {
				LogFile.writeLog("ERROR: agentTimeStep in agentGrid markup MUST be "+
						"less than or equal to the global timestep\n"+
						"\tagentTimeStep was given as: "+AGENTTIMESTEP+"\n"+
						"\tglobal time step is currently: "+SimTimer.getCurrentTimeStep());
				throw new Exception("agentTimeStep too large");
			}
			LogFile.writeLog("Agent time step is... " + value);
		}

		// Reference to the domain where this container is defined
		domain = (Domain) aSimulator.world.getDomain(root
				.getParam("computationDomain"));
		mySim = aSimulator;
		agentList = new LinkedList<SpecialisedAgent>();
		agentIter = agentList.listIterator();
		// Optimised the resolution of the grid used to sort located agents
		System.out.println("checkgridsize");
		checkGridSize(aSimulator, root);

		// Now initialise the padded grid

		createShovGrid(aSimulator);
		// Initialise spatial grid used to display species distribution
		createOutputGrid(aSimulator);

		if (Simulator.isChemostat){
			LogFile.writeLog("Chemostat volume is "+ExtraMath.cube(_res)+" micrometers cubed");
		}else{
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
	 * Call each grid cell of the agentGrid
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

			//sonia 27.04.2010
			//removing mislocalized/death agents from the agentList, so that we can proceed 
			//with the shoving algorithm
			// the piece of code below was taken from removeAllDead() but with the difference that we are not
			// removing the agents from the _agentToKill list (so that it can be used to record the
			//information on death agents later on)

			for(SpecialisedAgent aDeathAgent: _agentToKill){
				if (aDeathAgent.isDead) {
					nDead++;
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



		//SimTimer.setCurrentTimeStep(globalTimeStep);

		LogFile.chronoMessageOut("Agents stepped/dead/born: " + nAgent0 + "/"
				+ nDead + "/" + nBirth);

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
				try {
					removeOnBorder(this);
				} catch (Exception e) {
					System.out.println("At AgentContainer:removeOnBorder error met :" + e);
					System.exit(-1);
				}
				// mark biomass connected to the carrier and remove any non-connected portions
				if (DOSLOUGHING) {
					refreshGroupStatus();
					markForSloughing();
				}

				LogFile.chronoMessageOut("Detachment");
				
		}
		nAgent = agentList.size();
		if (maxPopLimit > 0 && nAgent >= maxPopLimit)
			aSim.continueRunning = false;
	}


	/**
	 * Compute pressure field and apply resulting advection movement
	 */
	public double followPressure() {
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
	 * Solve spatial spreading (acts only on located agents)
	 * 
	 * @param fullRelax
	 * @param maxShoveIter
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
	 * Used during initialisation to start from a coherent state
	 */
	public void relaxGrid() {

		//sonia:chemostat
		if(!Simulator.isChemostat){
			Collections.shuffle(agentList, ExtraMath.random);
			shoveAllLocated(true, true, MAXITER / 2, 0.1, 0.25);
			shoveAllLocated(true, true, MAXITER / 2, 0.1, 1);
		}
	}

	/* _________________________________________________________________ */
	/**
	 * 
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

	/* ____________________________ SHOVING FUNCTIONS _____________________ */

	protected void refreshGroupStatus() {
		for (int index = 0; index < _nTotal; index++)
			_grid[index].refreshElement();
	}

	/**
	 * Explore grid cells around the current one and sends (sonia: don't you mean returns?..) all agents contained :
	 * including grid cells on the other side of the cyclic boundary
	 * 
	 * @param agentPosition
	 * @param range :
	 *            maximal range to screen
	 * @param nbList:
	 *            the container where to store found particles
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

	public void registerBirth(SpecialisedAgent anAgent) {
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

	public void registerDeath(SpecialisedAgent anAgent) {

		//sonia 27.04.2010
		// just to make sure we are not adding death agents to the list that have already been added...
		if(!_agentToKill.contains(anAgent))
			_agentToKill.add(anAgent);

	}

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
	 * @author Sonia
	 *
	 * @param agentTimestep - this should be the same as the global timeStep or lower
	 * This method removes a number of agents from the system according to the dilution
	 * set for the chemostat. 
	 * The global time step should be set to 0.10*(1/D) so that around 10% of the agents
	 * will be removed from the system in each iteration.
	 * Remember, agents stand for all type of particles that can be removed (aka deleted) 
	 * from the system, from bacteria to eps.
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
			// Rob 13/09/2012: this method is now changed to include a running tally of the non-integer
			// remainder (as in removeOnBorder in biofilms)
			double temp = Dfactor*agentTimeStep*agentList.size() + tallyVariable;
			agentsToDilute = (int) Math.floor(temp);
			tallyVariable = temp%1;
			
		} else{
			// Rob 28/11/2011: Added this so we can keep the population at or below 1000
			// EROSIONMETHOD is normally used to decide between biofilm functions so this shouldn't be a problem
		    	if (Dfactor>0)
		    	    agentsToDilute = Math.max(agentList.size() - (int) Dfactor , 0);
		    	else
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
	 * 
	 * @param anAgent
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
	 * Update the position of the agent on the agent grid
	 * @param anAgent
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

	public void fitAgentMassOnGrid(SpatialGrid biomassGrid) {
		for (int i = 0; i < _nTotal; i++) {
			for (LocatedAgent aLoc : _grid[i].group) {
				aLoc.fitMassOnGrid(biomassGrid);
			}
		}
	}

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
	 * 
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
		for (int iSpecies = 0; iSpecies < aSim.speciesList.size(); iSpecies++) {
			_speciesGrid[iSpecies].writeReport(bufferState, bufferSum);
		}

		//		// output pressure field (not needed here b/c it's output w/ the solutes)
		//		if (aSim.getSolver("pressure").isActive())
		//			_pressure.writeReport(bufferState, bufferSum);

		/* Build a grid of detachment */
		//printLevelSet(bufferState);
	}

	/**
	 * 
	 * @param aSim
	 * @param bufferState
	 * @param bufferSum 
	 * @throws Exception
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
	 * @author Sonia
	 * recording data from agents that will be removed from the simulated environment
	 */
	public void writeReportDeath(Simulator aSim, ResultFile bufferStateDeath, ResultFile bufferSumDeath) throws Exception {


		if (!Idynomics.quietMode)
			System.out.println("size of agentToKill list at beginning of the writeReportDeath:  " + _agentToKill.size());

		// bvm 10.2.2009: include information about the shoving grid
		StringBuffer gridInfo = new StringBuffer();
		gridInfo.append("<grid");
		gridInfo.append(" resolution=\"").append(_res).append("\"");
		gridInfo.append(" nI=\"").append(_nI).append("\"");
		gridInfo.append(" nJ=\"").append(_nJ).append("\"");
		gridInfo.append(" nK=\"").append(_nK).append("\"");
		gridInfo.append("/>\n");
		//sonia 26.04.2010 writing result files for death/removed biomass
		bufferSumDeath.write(gridInfo.toString());
		bufferStateDeath.write(gridInfo.toString());

		//sonia: 5-05-09
		//include the information about the environment status at the beginning of the agent_State and
		//agent_Sum files
		if(Simulator.isFluctEnv){
			StringBuffer envInfo = new StringBuffer();
			envInfo.append("<Environment");
			envInfo.append(" env=\"").append(FluctEnv.envStatus).append("\"");
			envInfo.append("/>\n");
			bufferSumDeath.write(envInfo.toString());
			bufferStateDeath.write(envInfo.toString());
		}


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
					aSpecies.getProgenitor().sendHeader());
			//sonia 27.04.2010
			//introducing reason of death - header
			speciesBuffer[iSpecies].append(",death");
			speciesBuffer[iSpecies].append("\"");
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
		for (SpecialisedAgent anAgent : _agentToKill) {
			spIndex = anAgent.getSpecies().speciesIndex;
			spPop[spIndex]++;

			if (anAgent instanceof LocatedAgent) {
				aLoc = (LocatedAgent) anAgent;	
				spMass[spIndex] += aLoc.getTotalMass();
				spGrowth[spIndex] += aLoc.getNetGrowth();	

				//sonia:07-05-09
				/* <----- HGT Sonia Stats Start ----> */
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
				/* <----- HGT Sonia Stats End ----> */

				speciesBuffer[spIndex].append(aLoc.writeOutput());
				//sonia 27.04.2010
				//added death information to agent's state description
				speciesBuffer[spIndex].append("," + aLoc.death + ";\n");

			}
		}

		StringBuffer text;
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
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

			text.append("</species>\n");
			bufferSumDeath.write(text.toString());
		}

		//brian
		for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			speciesBuffer[iSpecies].append("</species>\n");
			bufferStateDeath.write(speciesBuffer[iSpecies].toString());
		}	


	}



	public void preprintLevelSet() {
		// Build the matrix of erosion time
		// _levelset.refreshFromBorder();

		for (int index = 0; index < _nTotal; index++) {
			_grid[index].printLevelSet(_erosionGrid);
		}
	}


	public void printLevelSet(ResultFile bufferState) throws Exception {

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
	}

	/* ________________ INITIALISATION __________________________ */

	public void checkGridSize(Simulator aSimulator, XMLParser root) {
		_res = root.getParamDbl("resolution");

		//sonia:chemostat
		if(Simulator.isChemostat){
			//set the resolution to the resolution of the domain
			_res = domain._resolution;
			//do not correct the grid size

		}else{
			// Eventually correct grid size
			_nI = (int) Math.ceil(domain.length_X / _res);
			_res = domain.length_X / _nI;

			_nJ = (int) Math.ceil(domain.length_Y / _res);
		}

		if (domain.is3D()) {
			_nK = (int) Math.ceil(domain.length_Z / _res);
			is3D = true;
		} else {
			_nK = 1;
			is3D = false;
		}

		//sonia:chemostat
		if (Simulator.isChemostat){
			//sonia:chemostat
			//set the number of grid elements to 1
			_nTotal = 1;
		}else{
			_nTotal = (_nI + 2) * (_nJ + 2) * (_nK + 2);
		}

	}

	/**
	 * Create a vectorized array of spatial groups, build their neighbourhood
	 * and store it Note :Shoving grid is a padded of LocatedGroups
	 */
	public void createShovGrid(Simulator aSimulator) {

		_grid = new LocatedGroup[_nTotal];	
		for (int index = 0; index < _nTotal; index++)
			_grid[index] = new LocatedGroup(index, this, aSimulator);

		for (int index = 0; index < _nTotal; index++){
			_grid[index].init();}


	}

	public void createOutputGrid(Simulator aSim) {
		_erosionGrid = new double[_nI + 2][_nJ + 2][_nK + 2];
		_speciesGrid = new SpatialGrid[aSim.speciesList.size()];
		for (int iSpecies = 0; iSpecies < aSim.speciesDic.size(); iSpecies++)
			_speciesGrid[iSpecies]=domain.createGrid(aSim.speciesDic.get(iSpecies), 0);
	}

	/* _____________________ EROSION & DETACHMENT ___________________________ */

	/**
	 * 
	 */
	protected void markForSloughing() {
		// perform connected volume filtration (connected to bottom
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
	 * 
	 * 
	 */
	public void shrinkOnBorder() {
		// find the border points (erosion and sloughing have changed the
		// configuration)


		_levelset.refreshBorder(true, mySim);

		double mass = 0;
		double tallyVariable = 0;
		int nDetach = 0;
		int index;
		double ratio;

		// System.out.print("erosion ratio: ");
		for (LocatedGroup aBorderElement : _levelset.getBorder()) {
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

		//sonia 26.04.2010
		//commented out removeAllDead();
		//removeAllDead();

		LogFile.writeLog("Eroding " + nDetach + " ("
				+ ExtraMath.toString(mass, true) + "/"
				+ ExtraMath.toString(tallyVariable, true) + " fg)");
	}
	

	private double detFunction(double i, double location) {
		if (i<1){ return ExtraMath.sq(_res - (location % _res));}
		else {return ExtraMath.sq(location % _res);}
	}

	/**
	 * This is a variant on OLDremoveOnBorder(), which performs the discrete removal of agents over the whole 
	 * _close group together, instead of considering each Located Group seperately. Note also that here tallyVariable
	 * is cumulative: whatever remains rolls over to the next time-step. This eliminates any timestep issues.
	 * @param agentGrid
	 */
	public void removeOnBorder(AgentContainer agentGrid) {
	// public void NEWremoveOnBorder(AgentContainer agentGrid) {
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

		LogFile.writeLog("Eroding " + nDetach + " ("
				+ ExtraMath.toString(mass, true) + "/"
				+ ExtraMath.toString(tallyVariable, true) + " fg) from "
				+ _levelset.getBorder().size() +" elements.");
	}


	public void calcDetPriority(AgentContainer agentGrid, LocatedGroup aBorderElement, double ratio){
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
	 * Check that these coordinates are defined on this grid (the grid is padded
	 * but dc uses shifted coordinates)
	 */
	public boolean isValid(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i >= 0) & (dC.i < _nI);
		out &= (dC.j >= 0) & (dC.j < _nJ);
		out &= (dC.k >= 0) & (dC.k < _nK);
		return out;
	}

	// Check that the nb grid cell has valid coordinates
	public boolean isValid(ContinuousVector cC) {
		return isValid(getGridPosition(cC));
	}



	/**
	 * find the voxel a continuous coordinate lies in and return the index
	 * 
	 * @param cc:the
	 *            continuous coordinate to find the index for
	 * @return index on the 1D (vectorized) array
	 */
	public int getIndexedPosition(ContinuousVector cc) {

		//sonia:chemostat
		// to guarantee that the agents are effectively removed
		if(Simulator.isChemostat){
			return 0;

		}else{

			int i = (int) Math.floor(cc.x / _res) + 1;
			int j = (int) Math.floor(cc.y / _res) + 1;
			int k = (int) Math.floor(cc.z / _res) + 1;

			return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
		}
	}


	public int getIndexedPosition(DiscreteVector dc) {
		int i = dc.i + 1;
		int j = dc.j + 1;
		int k = dc.k + 1;

		return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
	}

	public DiscreteVector getGridPosition(ContinuousVector cC) {
		int i = (int) Math.floor(cC.x / _res);
		int j = (int) Math.floor(cC.y / _res);
		int k = (int) Math.floor(cC.z / _res);

		return new DiscreteVector(i, j, k);
	}

	public DiscreteVector getGridPosition(int index) {
		int k = (int) Math.floor(index / (_nI + 2) / (_nJ + 2));
		int j = (int) Math.floor((index - k * ((_nI + 2) * (_nJ + 2)))
				/ (_nI + 2));
		int i = (int) Math.floor((index - (k * (_nI + 2) * (_nJ + 2)) - j
				* (_nI + 2)));

		return new DiscreteVector(i - 1, j - 1, k - 1);
	}

	public ContinuousVector getGridLocation(int index) {
		int k = (int) Math.floor(index / (_nI + 2) / (_nJ + 2));
		int j = (int) Math.floor((index - k * ((_nI + 2) * (_nJ + 2)))
				/ (_nI + 2));
		int i = (int) Math.floor((index - (k * (_nI + 2) * (_nJ + 2)) - j
				* (_nI + 2)));

		return new ContinuousVector((i + .5 - 1) * _res, (j + .5 - 1) * _res,
				(k + .5 - 1) * _res);
	}



	public double getAgentTimeStep() {
		return AGENTTIMESTEP;
	}


	public double getResolution() {
		return _res;
	}

	public int[] getGridDescription() {
		int[] out = { _nI, _nJ, _nK };
		return out;
	}

	public LocatedGroup[] getShovingGrid() {
		return _grid;
	}

	public LevelSet getLevelSet() {
		return _levelset;
	}

}
