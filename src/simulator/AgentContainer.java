/**
 * \package simulator
 * \brief Package of classes that create a simulator object and capture
 * simulation time.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator;

import java.util.*;

import idyno.SimTimer;
import simulator.agent.*;
import simulator.agent.zoo.MultiEpiBac;
import simulator.agent.zoo.MultiEpisome;
import simulator.detachment.*;
import simulator.diffusionSolver.DiffusionSolver;
import simulator.diffusionSolver.Solver_pressure;
import simulator.geometry.*;
import simulator.geometry.boundaryConditions.ConnectedBoundary;
import simulator.SpatialGrid;
import utils.ResultFile;
import utils.XMLParser;
import utils.LogFile;
import utils.ExtraMath;

/**
 * \brief Class to store all the agents, call them, and manage shoving/erosion
 * of located agents.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public class AgentContainer 
{
	/**
	 * Computational domain to which this grid is assigned
	 */ 
	public Domain domain;
	
	/**
	 * Local copy of the simulation object used to create the conditions
	 * specified in the protocol file.
	 */
	public Simulator mySim;

	/**
	 * Container for all agents (even the non located ones)
	 */
	public LinkedList<SpecialisedAgent> agentList;
	
	/**
	 * Temporary containers used to store agents who will be added or removed.
	 * Visibility public so that it can be accessed from LocatedGroup in
	 * killAll().
	 */
	public LinkedList<SpecialisedAgent> _agentToKill = 
										new LinkedList<SpecialisedAgent>();

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
	 * Resolution of the grid. Specified in the XML protocol file.
	 */
	private Double _res;
	
	
	/**
	 *	Grid to hold the agents in this container. This grid holds groups of
	 *agents, within a LocatedGroup object. 
	 */
	private LocatedGroup[] _grid;
	
	/**
	 * 3D array that captures erosion in the agent grid
	 */
	protected Double[][][] _erosionGrid;
	
	/**
	 * Boolean noting whether the grid is 3D or 2D
	 */ 
	public boolean is3D;
	
	/**
	 * Sets the allowed cut off for shoving agents. When less that this
	 * fraction of agents are still moving due to shoving, stop. 
	 */
	private final Double SHOVEFRACTION;

	/**
	 * Allows you to define a smaller timestep for agent behaviors and
	 * interactions. Read in by Simulator class. This is a local copy. 
	 */
	public final Double AGENTTIMESTEP;
	
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
	double tallyVariable = 0.0; 
	
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
		SHOVEFRACTION = root.getParamDbl("shovingFraction");
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

		double value = agentTimeStep;

		// Now deal with the agent timestep
		if (Double.isNaN(value))
		{
			AGENTTIMESTEP = SimTimer.getCurrentTimeStep();
			LogFile.writeLog("Using global timestep of "+AGENTTIMESTEP+" for agentTimeStep");
		}
		else
		{
			AGENTTIMESTEP = value;
			if (AGENTTIMESTEP > SimTimer.getCurrentTimeStep()) 
			{
				LogFile.writeLog("ERROR: agentTimeStep in agentGrid markup MUST be "+
						"less than or equal to the global timestep\n"+
						"\tagentTimeStep was given as: "+AGENTTIMESTEP+"\n"+
						"\tglobal time step is currently: "+SimTimer.getCurrentTimeStep());
				throw new Exception("agentTimeStep too large");
			}
			LogFile.writeLog("Agent time step is... " + value);
		}

		// Now set the domain where this container is defined
		domain = (Domain) aSimulator.world.getDomain(root.getParam("computationDomain"));
		mySim = aSimulator;
		
		agentList = new LinkedList<SpecialisedAgent>();
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
			_levelset = LevelSet.staticBuilder(root.getChildParser("detachment"), this);

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
	public void step(Simulator aSim)
	{
		SpecialisedAgent anAgent;
		/* STEP AGENTS ________________________________________________ */
		LogFile.chronoMessageIn();
		Collections.shuffle(agentList, ExtraMath.random);
		
		// Record values at the beginning
		int nBirth = 0;
		int nAgent = agentList.size();
		double dt = 0.0;
		double elapsedTime = 0.0;
		double globalTimeStep = SimTimer.getCurrentTimeStep();
		// for the local time step, choose the value according to which is best
		double localdt = Math.min(AGENTTIMESTEP,globalTimeStep);

		int nAgent0 = agentList.size();
		// Apply a shorter time step when visiting all the agents

		while (elapsedTime < globalTimeStep)
		{
			// by default use the saved agent timestep
			dt = localdt;


			// check for a smaller dt (usually for the last iterate)
			if (dt > (globalTimeStep-elapsedTime))
				dt = globalTimeStep-elapsedTime;

			elapsedTime += dt;		

			/* Step all the agents */
			SimTimer.setCurrentTimeStep(dt);

			// Bypass agent movement in a chemostat.
			if ( ! Simulator.isChemostat )
				followPressure();
			
			for ( SpecialisedAgent agent : agentList )
				agent.step();
			
			Collections.shuffle(agentList, ExtraMath.random);

			if ( Simulator.isChemostat )
				agentFlushedAway(dt);


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
			
			// Apply moderate overlap relaxation, unless this is a chemostat.
			if( ! Simulator.isChemostat )
				shoveAllLocated(15);
			
		}

		SimTimer.setCurrentTimeStep(globalTimeStep);
		
		
		// KA - MOVED OUTPUT OF AGENTS STEPPED / DEAD / BORN FROM HERE TO LATER, SO THAT WE CAN INCLUDE ERODED CELLS IN THE
		// COUNT OF DEAD CELLS
				
		/* MECHANICAL INTERACTIONS _____________________________________ */
		
		if( ! Simulator.isChemostat )
		{
			//sonia 26.04.2010
			//care as been take so that the death agents are removed from the 
			// _agentList preventing their participation in the shoving process 			
			// spring and then shove only particles
			LogFile.chronoMessageIn("Shoving");
			shoveAllLocated(MAXITER);
			LogFile.chronoMessageOut("Shoving done");
			
			
			// EROSION & DETACHMENT _________________________________________ */
			// Refresh the space occupation map (-1:outside, 0:carrier,1:biofilm, 2:liquid, 3:bulk)
			LogFile.chronoMessageIn("Detachment/erosion");
			refreshGroupStatus();
			/*
			 * Rebuild the border of the biofilm and compute erosion-time for
			 * the whole biofilm
			 */
			_levelset.refreshBorder(true, mySim);
			_levelset.computeLevelSet(mySim);
			/*
			 * On grid elements on the border apply a probabilistic erosion
			 */
			// Rob Feb 2011: added borderErosion(), an alternative to erodeBorder(), which
			// removes whole agents in a deterministic way
			if (EROSIONMETHOD)
				shrinkOnBorder();
			else
				// KA - Check added after self-attachment, as it is possible that if the timestep is much less than the input rate, 
				// the grid may contain no cells for the first few steps
				if ( ! this.agentList.isEmpty() )
				{
					try
					{
						LogFile.writeLogAlways("Trying to removeOnBorder");
						removeOnBorder(this);
					}
					catch (Exception e)
					{
						LogFile.writeError(e, "AgentContainer.removeOnBorder()");
						System.exit(-1);
					}
					// mark biomass connected to the carrier and remove any non-connected portions
					if (DOSLOUGHING) {
						refreshGroupStatus();
						markForSloughing();
					}
				}

			LogFile.chronoMessageOut("Detachment/erosion done");
			
				
		}
		
		// OUTPUT THE COUNT STATISTICS
		LogFile.chronoMessageOut("Agents stepped/dead/born: " + nAgent0 + "/"
				+ _agentToKill.size() + "/" + nBirth);

		
		nAgent = agentList.size();
		if (maxPopLimit > 0 && nAgent >= maxPopLimit)
			aSim.continueRunning = false;
	}


	/**
	 * \brief Compute pressure field and apply resulting advection movement to
	 * affected agents.
	 */
	public void followPressure() 
	{
		DiffusionSolver solver = mySim.getSolver("pressure");
		
		// Find a solver for pressure field and use it
		// don't use the pressure if it's not active
		if ( solver == null || ! solver.isActive() )
			return;
		
		LogFile.writeLog("Doing pressure calculations.");
		
		// get local timestep (which was set in the step() routine calling this one)
		Double dt = SimTimer.getCurrentTimeStep();
		
		// Solve for pressure field
		solver.initAndSolve();
		_pressure = ((Solver_pressure) solver).getPressureGrid();
		
		// copy calculated pressure field to the solute list
		// (allows easy output of pressure field)
		mySim.getSolute("pressure").setGrid(_pressure.getGrid());

		// Determine local advection speed
		Double maxSpeed = 0.0;
		for (LocatedGroup aGroup : _grid)
		{
			maxSpeed = Math.max(maxSpeed,
								aGroup.computeMove(_pressure, AGENTTIMESTEP));
		}
		
		// bvm 04.03.09: new method to address any high velocities:
		// use smaller local timesteps to keep the movement under control
		Double dtlocal = dt;
		int itlocal = 1;
		
		while ( maxSpeed > this._res/dtlocal )
		{
			// if the move takes an agent farther than one grid element,
			// apply scaling factor until move is within limit
			dtlocal /= 10.0;
			itlocal *= 10;
		}
		
		if (itlocal > 1)
		{
			LogFile.writeLog("PRESSURE MOVEMENT HAS LOCAL TIMESTEP "
					+dtlocal+" ("+itlocal+" iterations)");
		}
		
		// scale movement vectors based on new, smaller timestep and apply
		// the movement to each agent in each group
		Double alpha = dtlocal/dt;
		for (LocatedGroup aGroup : _grid)
			aGroup.addMoveToAgents(alpha);
		
		// now apply the scaled agent movements to each agent
		for (int i = 0; i < itlocal; ++i)
			for ( SpecialisedAgent anAgent : agentList )
				anAgent.move();
	}	


	/**
	 * \brief Solve spatial spreading through application of shoving (acts
	 * only on located agents).
	 * 
	 * @param maxShoveIter	The maximum number of shoving iterations that
	 * should be applied to find a new position.
	 */
	public void shoveAllLocated(int maxShoveIter)
	{
		int nMoved;
		shovLimit = Math.max(1, (int) (agentList.size() * SHOVEFRACTION));
		shovIter = 0;
		do 
		{
			nMoved = performMove();
		} while ((shovIter++ < maxShoveIter) && (nMoved >= shovLimit));
		LogFile.writeLog(nMoved + "/" + agentList.size() + " after " + shovIter
				+ " shove iterations");
	}

	/**
	 * \brief Used during initialization to start from a coherent state.
	 * 
	 * TODO Rob 13Mar2015: Why are we doing this five times the "maximum"?
	 */
	public void relaxGrid() 
	{
		if( ! Simulator.isChemostat )
		{
			Collections.shuffle(agentList, ExtraMath.random);
			shoveAllLocated(5 * MAXITER);
		}
	}


	/**
	 * \brief Moves an agent as a result of shoving or pushing due to growth
	 *
	 * @param isSynchro
	 */
	protected int performMove()
	{
		int nMoved = 0;
		Double deltaMove;
		/*
		 * Compute movement, deltaMove is relative movement.
		 */
		for ( SpecialisedAgent agent : agentList )
		{
			deltaMove = agent.interact(MUTUAL);
			nMoved += (deltaMove >= 0.1  ? 1 : 0);
		}
		return nMoved;
	}


	/**
	 * \brief Refresh the space occupation map as agents may have moved
	 * 
	 * Each grid square has a status:
	 * -1 outside
	 * 0  carrier
	 * 1  biofilm 
	 * 2  liquid 
	 * 3  bulk
	 */
	protected void refreshGroupStatus()
	{
		for ( LocatedGroup aLG : _grid )
			aLG.refreshElement();
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
	public void getPotentialShovers(int index, Double range,
											LinkedList<LocatedAgent> nbList)
	{
		LocatedGroup aGroup;
		int radius = Math.max(1, (int) Math.floor(range / this._res));
		nbList.clear();
		for (int i = -radius; i <= radius; i++)
		{
			if ( _grid[index].moveX(i) == null )
				continue;
			for (int j = -radius; j <= radius; j++)
			{
				if ( is3D )
				{
					if (_grid[index].moveX(i).moveY(j) == null)
						continue;
					for (int k = -radius; k <= radius; k++)
					{
						aGroup = _grid[index].moveX(i).moveY(j).moveZ(k);
						if ( aGroup != null )
							nbList.addAll(aGroup.group);
					}
				}
				else
				{
					aGroup = _grid[index].moveX(i).moveY(j);
					if ( aGroup != null )
						nbList.addAll(aGroup.group);
				}
			}
		}
	}

	/* ________________ TOOLS:GRID, MAP & TREE MANAGEMENT __________________ */

	/**
	 * \brief Registers the birth of an agent and adds this to the agent grid.
	 * 
	 * @param anAgent	New agent to add to the agent grid.
	 */
	public void registerBirth(SpecialisedAgent anAgent) 
	{
		// Add the agent to agentList
		agentList.add(anAgent);

		// Add the agent on the grid
		if (anAgent instanceof LocatedAgent)
		{
			LocatedAgent aLoc = (LocatedAgent) anAgent;
			try
			{
				if ( Simulator.isChemostat )
					_grid[0].add(aLoc);
				else
					_grid[getIndexedPosition(aLoc.getLocation())].add(aLoc);
			}
			catch (Exception e)
			{
				LogFile.writeError(e, "AgentContainer.registerBirth()");
			}
		}

	}

	/**
	 * \brief Register the death of an agent, by adding it to a list of agents
	 * that will be removed from the simulation.
	 * 
	 * @param anAgent	SpecialisedAgent object that will be removed from the
	 * simulation.
	 */
	public void registerDeath(SpecialisedAgent anAgent) 
	{
		if ( ! _agentToKill.contains(anAgent) )
			_agentToKill.add(anAgent);
	}

	/**
	 * \brief Iterates through the _agentToKill list and removes these agents
	 * from the grid. 
	 * 
	 * These cells should be assumed to be dead. Returns the number of dead
	 * cells removed.
	 * 
	 * @return	Number of dead cells removed from the simulation.
	 */
	public int removeAllDead()
	{
		int nDead = 0;
		ListIterator<SpecialisedAgent> iter = _agentToKill.listIterator();
		SpecialisedAgent anAgent;

		while (iter.hasNext())
		{
			anAgent = iter.next();
			if (anAgent.isDead)
			{
				nDead++;
				iter.remove();
				agentList.remove(anAgent);
				removeLocated(anAgent);
			}
		}
		
		/* TODO Rob 13Mar2015: Started tidying up this.
		for ( SpecialisedAgent aSA : _agentToKill )
		{
			if ( ! aSA.isDead )
			{
				_agentToKill.remove(aSA);
				continue;
			}
			removeLocated(aSA);
		}
		
		nDead = _agentToKill.size();
		agentList.removeAll(_agentToKill);
		*/
		
		_agentToKill.clear();
		return nDead;
	}

	/**
	 * \brief Removes a number of agents from the system according to the
	 * dilution set for the chemostat.
	 * 
	 * The global time step should be set to 0.10*(1/D) so that around 10% of
	 * the agents will be removed from the system in each iteration. Remember,
	 * agents stand for all type of particles that can be removed (i.e.
	 * deleted) from the system, from bacteria to EPS.
	 *
	 * @param agentTimeStep - this should be the same as the global timeStep
	 * or lower.
	 */
	public void agentFlushedAway(Double agentTimeStep)
	{
		/*
		 * After having shuffled the list (during the step()) with all the
		 * agents we are now ready to kill agents according to the dilution
		 * value read from the Bulk class.
		 * 
		 * Sonia:2.03.2010 Just in case, let's do it again
		 */
		Collections.shuffle(agentList, ExtraMath.random);
		Bulk aBulk;
		for ( ConnectedBoundary aBC : domain.getAllConnectedBoundaries() )
		{
			aBulk = ((ConnectedBoundary) aBC).getBulk();
			if ( aBulk.nameEquals("chemostat") )
				Dfactor = aBulk._D;
		}
		
		int agentsToDilute = 0;
		if (EROSIONMETHOD)
		{
			Double temp = Dfactor*agentTimeStep*agentList.size() + tallyVariable;
			agentsToDilute = temp.intValue();
			tallyVariable = temp % 1;
		}
		else
		{
			/*
			 * Rob 28/11/2011: Added this so we can keep the population at or
			 * below 1000. EROSIONMETHOD is normally used to decide between
			 * biofilm functions so this shouldn't be a problem.
			 * 
			 * TODO Rob 16Mar2015: Make more robust.
			 */
			agentsToDilute = Math.max(agentList.size() - 1000, 0);
		}
		
		for (SpecialisedAgent anAgent : agentList.subList(0, agentsToDilute))
		{
			anAgent.isDead = true;
			anAgent.death = "dilution";
			anAgent.die(false);
		}
		
		// TODO Rob 13Mar2015: Simplify? agentList.removeAll(_agentToKill);
		for ( SpecialisedAgent anAgent : agentList )
			if ( anAgent.isDead )
				agentList.remove(anAgent);
	}
	
	/**
	 * \brief Remove an agent from the grid.
	 * 
	 * @param anAgent	Agent to be removed from the grid.
	 */
	public void removeLocated(SpecialisedAgent anAgent)
	{
		if (anAgent instanceof LocatedAgent)
		{
			LocatedAgent aLoc = (LocatedAgent) anAgent;
			int index = getIndexedPosition(aLoc.getLocation());
			if ( ! Double.isNaN(index) )
				_grid[index].remove(aLoc);
		}
	}
	
	/**
	 * \brief Update the position of the agent on the agent grid.
	 * 
	 * The agent may have moved due to shoving, pushing, or growth.
	 * 
	 * @param anAgent	Agent that has moved and needs the location updated.
	 */
	public void registerMove(LocatedAgent anAgent)
	{
		/*
		 * Compute the theoretical index on the agentGrid
		 */
		int newIndex = getIndexedPosition(anAgent.getLocation());
		int oldIndex = anAgent.getGridIndex();
		/*
		 * If gridIndex has changed, update the references.
		 */
		LogFile.writeLogDebug("Debugging AgentContainer.registerMove()");
		if ( isValid(anAgent.getLocation()) )
		{
			LogFile.writeLogDebug("Agent location "+
								anAgent.getLocation().toString()+"is valid");
			if ( newIndex != oldIndex )
			{
				_grid[oldIndex].remove(anAgent);
				_grid[newIndex].add(anAgent);
				anAgent.setGridIndex(newIndex);
			}
		}
		else
		{
			LogFile.writeLogDebug("Agent location "+
				anAgent.getLocation().toString()+" is not valid -> Killed");
			//anAgent.death = "overBoard";
			anAgent.die(false);
		}
	}
	
	/**
	 * \brief Iterates through each square of the agent grid to calculate
	 * total biomass in that grid space, storing this on the supplied biomass
	 * grid.
	 * 
	 * @param biomassGrid	The biomass grid that will contain the total
	 * biomass in each grid square of the agent grid.
	 */
	public void fitAgentMassOnGrid(SpatialGrid biomassGrid) 
	{
		for ( LocatedGroup aSquare : _grid )
			for ( LocatedAgent aLoc : aSquare.group )
				aLoc.fitMassOnGrid(biomassGrid);
	}

	/**
	 * \brief Iterates through each square of the agent grid to calculate
	 * total volume in that grid space, storing this on the supplied biomass
	 * grid.
	 * 
	 * @param biomassGrid	The biomass grid that will contain the total
	 * volume in each grid square of the agent grid.
	 */
	public void fitAgentVolumeRateOnGrid(SpatialGrid biomassGrid)
	{
		biomassGrid.resetToZero();
		for (LocatedGroup aSquare : _grid)
			for (LocatedAgent aLoc : aSquare.group)
				aLoc.fitVolRateOnGrid(biomassGrid);
	}



	/* ____________________ REPORT FILE EDITION__________________________ */

	/**
	 * \brief Summarize the agents in this container in simulation output, the
	 * form of a total biomass grid.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param bufferState	The agent_state result file output buffer.
	 * @param bufferSum	The agent_sum result file output buffer.
	 * @throws Exception	Exception thrown if there are issues writing to
	 * these buffers.
	 */
	public void writeGrids(Simulator aSim, ResultFile bufferState,
									ResultFile bufferSum) throws Exception
	{
		LocatedAgent aLoc;

		//sonia:chemostat
		//I've modified the refreshElement() method for the chemostat case

		// Refresh the space occupation map (0:carrier,1:biofilm or 2:liquid)
		for (LocatedGroup aGroup : _grid)
			aGroup.refreshElement();
		
		/* Build a grid of biomass concentration */
		
		// Set existing grid to zero
		for (SpatialGrid aSpeciesGrid : _speciesGrid)
			aSpeciesGrid.resetToZero();
		
		// Sum biomass concentrations
		for (SpecialisedAgent anA : agentList)
			if (anA instanceof LocatedAgent)
			{
				aLoc = (LocatedAgent) anA;
				aLoc.fitMassOnGrid(_speciesGrid[aLoc.speciesIndex]);
			}

		// now output the biomass values
		for (SpatialGrid aSpeciesGrid : _speciesGrid)
			aSpeciesGrid.writeReport(bufferState, bufferSum);
	}

	/**
	 * \brief Write simulation output to XML files, summarizing information
	 * for each species and plasmid on the agent grid.
	 * 
	 * These files are written out as
	 * agent_state([simulation step number]).xml and
	 * agent_sum([simulation step number]).xml
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param bufferState	The agent_state result file output buffer.
	 * @param bufferSum	The agent_sum result file output buffer.
	 * @throws Exception	Exception thrown if there are issues writing to
	 * these buffers.
	 */
	public void writeReport(Simulator aSim, ResultFile bufferState,
										ResultFile bufferSum) throws Exception
	{
		/*
		 * This will be our general-purpose buffer.
		 */
		StringBuffer textBuffer = new StringBuffer();
		/*
		 * Get the number of species in this simulation.
		 */
		int nSpecies = aSim.speciesList.size();
		/*
		 * Set up a buffer to hold information for each agent of these
		 * species.
		 */
		StringBuffer[] speciesBuffer = new StringBuffer[nSpecies];
		/*
		 * Include information about the shoving grid.
		 */
		textBuffer = writeGridInformation(textBuffer);
		/*
		 * Include the information about the environment status at the
		 * beginning of the agent_State and agent_Sum files.
		 */
		if ( Simulator.isFluctEnv )
			textBuffer = writeFluctEnvInformation(textBuffer);
		/*
		 * Write the header to the agent_State and agent_Sum buffers.
		 */
		bufferState.write(textBuffer);
		bufferSum.write(textBuffer);
		/*
		 * Initialise a Species markup for each present species
		 */
		for (Species aSpec : aSim.speciesList)
		{
			textBuffer  = new StringBuffer();
			/*
			 * First part, the species name.
			 */
			textBuffer.append("<species name=\"");
			textBuffer.append(aSpec.speciesName);
			/*
			 * Now create the header. Note that this information comes from
			 * the species itself - its part of the agent and the inherited
			 * classes. So if we are going to change that, it gets changed
			 * there.
			 */
			textBuffer.append("\" header=\"");
			textBuffer.append(aSpec.getProgenitor().sendHeader());
			textBuffer.append("\" >\n");
			speciesBuffer[aSpec.speciesIndex] = textBuffer;
		}
		/*
		 *  Initialise statistics (population total mass, growth-rate).
		 */
		int[] spPop = new int[nSpecies];
		Double[] spMass = ExtraMath.newDoubleArray(nSpecies);
		Double[] spGrowth = ExtraMath.newDoubleArray(nSpecies);
		
		/* <----- HGT Stats Begin ---------> */
		int plasmidListSize = 0;
		if ( Simulator.multiEpi )
			plasmidListSize = aSim.plasmidList.size();
		int[] spPlasmid = new int[nSpecies];
		int[] spConjugEvents = new int[nSpecies];
		int[][] spPlasmidTypes = new int [nSpecies][plasmidListSize];
		/* <----- HGT Stats End ----> */
		
		// Fill the agent_state file, build the state for the summary
  		LocatedAgent aLoc;
  		MultiEpiBac anEpiBac;
 		int spIndex;
 		for (SpecialisedAgent anAgent : agentList)
 		{
 			spIndex = anAgent.getSpecies().speciesIndex;
 			spPop[spIndex]++;
 			// Skip to the next agent if this one is dead
 			// TODO RC - do we really want to include dead agents in the
 			// population count?
 			if (anAgent.isDead)
 				continue;
 			
 			// TODO RC - why aren't we including the mass and growth rates
 			// of all ActiveAgents, only LocatedAgents?
 			if (anAgent instanceof LocatedAgent)
 			{
 				aLoc = (LocatedAgent) anAgent;	
 				spMass[spIndex] += aLoc.getTotalMass();
 				spGrowth[spIndex] += aLoc.getNetGrowth();
  				speciesBuffer[spIndex].append(aLoc.writeOutput()+";\n");
 			}
 			
 			/*<-------HGT Sonia Begin ------> */
 			// For each MultiEpiBac that currently hosts at least one plasmid:
 			if(anAgent instanceof MultiEpiBac )
 			{
 				anEpiBac = (MultiEpiBac) anAgent;
 				if (anEpiBac.plasmidHosted.size() > 0)
 				{
 					// Let 
 					spPlasmid[spIndex]++;
 					for (MultiEpisome epiPlasmid : anEpiBac.plasmidHosted)
 						for (String simPlasmid : aSim.plasmidList)
 							if (epiPlasmid.getName().equals(simPlasmid))
 								spPlasmidTypes[spIndex]
 									[aSim.plasmidList.indexOf(simPlasmid)]++;
 				}
 			 	if( anEpiBac.plasmidVector.size() > 0 )
 			 		spConjugEvents[spIndex] += anEpiBac.plasmidVector.size();
 			}
 			/* <------- HGT Sonia End ------> */
 			
 		}
 		
 		for (Species aSpecies : aSim.speciesList)
 		{
 			spIndex = aSpecies.speciesIndex;
 			// Write the agent_state info for this species to file  
 			speciesBuffer[spIndex].append("</species>\n");
 			bufferState.write(speciesBuffer[spIndex]);
 			
 			// Collate the agent_Sum info.
 			textBuffer = new StringBuffer();
 			textBuffer.append("<species name=\"");
 			textBuffer.append(aSpecies.speciesName);
 			textBuffer.append("\" header=\"population,mass,growthRate");
 			
 			
 			/*<----HGT Sonia begin---->*/
 			if (Simulator.multiEpi)
 			{
 				textBuffer.append(",plasmidBearing");
 				if (plasmidListSize > 0)
 					for (String plName : aSim.plasmidList)
 						textBuffer.append(",").append(plName);
 				textBuffer.append(",conjCount");
 			}
 			/*<----HGT Sonia end---->*/
 			
 			textBuffer.append("\" >\n");
 			
 			textBuffer.append(spPop[spIndex]).append(",");
 			textBuffer.append(spMass[spIndex]).append(",");
 			textBuffer.append(spGrowth[spIndex]);
 			
 			/*<----HGT Sonia begin---->*/
 			if(Simulator.multiEpi)
 			{
 				textBuffer.append(",").append(spPlasmid[spIndex]);
 				for (int c=0; c<spPlasmidTypes[spIndex].length; c++)
 					textBuffer.append(",").append(spPlasmidTypes[spIndex][c]);
 				textBuffer.append(",").append(spConjugEvents[spIndex]);
 			}
 			/*<----HGT Sonia end---->*/
 			
 			textBuffer.append(";\n</species>\n");
 			bufferSum.write(textBuffer);
 			
 		}
	}


	/**
	 * \brief Output information on agents that have been removed from the
	 * simulation.
	 * 
	 * Will enable investigators to find out why this is the case.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param bufferStateDeath	The agent_state_death result file output
	 * buffer.
	 * @param bufferSumDeath	The agent_sum_death result file output buffer.
	 * @throws Exception	Exception thrown if there are issues writing to
	 * these buffers.
	 */
	public void writeReportDeath(Simulator aSim, ResultFile bufferStateDeath,
								ResultFile bufferSumDeath) throws Exception
	{
		/*
		 * Status message to screen (not in quiet mode)
		 */
		LogFile.writeLog("size of agentToKill list at beginning of the"+
								"writeReportDeath: "+ _agentToKill.size());
		/*
		 * This will be our general-purpose buffer.
		 */
		StringBuffer textBuffer = new StringBuffer();
		/*
		 * Get the number of species in this simulation.
		 */
		int nSpecies = aSim.speciesList.size();
		/*
		 * Set up a buffer to hold information for each agent of these
		 * species.
		 */
		StringBuffer[] speciesBuffer = new StringBuffer[nSpecies];
		/*
		 *  Include information about the shoving grid.
		 */
		textBuffer = writeGridInformation(textBuffer);
		/*
		 * Include the information about the environment status at the
		 * beginning of the agent_StateDeath and agent_SumDeath files.
		 */
		if ( Simulator.isFluctEnv )
			textBuffer = writeFluctEnvInformation(textBuffer);
		/*
		 * Write the header to the agent_StateDeath and agent_SumDeath
		 * buffers.
		 */
		bufferStateDeath.write(textBuffer);
		bufferSumDeath.write(textBuffer);
		/*
		 * Initialise a Species markup for each present species
		 */
		for (Species aSpec : aSim.speciesList)
		{
			textBuffer  = new StringBuffer();
			/*
			 * First part, the species name.
			 */
			textBuffer.append("<species name=\"");
			textBuffer.append(aSpec.speciesName);
			/*
			 * Now create the header. Note that this information comes from
			 * the species itself - its part of the agent and the inherited
			 * classes. So if we are going to change that, it gets changed
			 * there.
			 */
			textBuffer.append("\" header=\"");
			textBuffer.append(aSpec.getProgenitor().sendHeader());
			/*
			 * Now we can add additional fields if these are not declared by
			 * the agent itself. In this case, we're going to add a reason
			 * for the agent's death.
			 */
			textBuffer.append(",death");
			textBuffer.append("\" >\n");
			speciesBuffer[aSpec.speciesIndex] = textBuffer;
		}
		/*
		 * Initialise statistics (population total mass, growth-rate).
		 */
		int[] spPop = new int[nSpecies];
		Double[] spMass = ExtraMath.newDoubleArray(nSpecies);
		Double[] spGrowth = ExtraMath.newDoubleArray(nSpecies);
		
		/* <----- HGT Stats Begin ---------> */
		int plasmidListSize = 0;
		if ( Simulator.multiEpi )
			plasmidListSize = aSim.plasmidList.size();
		int[] spPlasmid = new int[nSpecies];
		int[] spConjugEvents = new int[nSpecies];
		int[][] spPlasmidTypes = new int [nSpecies][plasmidListSize];
		/* <----- HGT Stats End ----> */
		/*
		 * Collate the information for the agent_StateDeath file.
		 */
		LocatedAgent aLoc;
  		MultiEpiBac anEpiBac;
  		int spIndex;
  		for (SpecialisedAgent anAgent : _agentToKill)
  		{
  		  	spIndex = anAgent.getSpecies().speciesIndex;
  		  	spPop[spIndex]++;
  		  	
  		  	// TODO RC - why aren't we including the mass and growth rates
  			// of all ActiveAgents, only LocatedAgents?
  			if (anAgent instanceof LocatedAgent)
  			{
  				aLoc = (LocatedAgent) anAgent;
  				spMass[spIndex] += aLoc.getTotalMass();
  				spGrowth[spIndex] += aLoc.getNetGrowth();
  				speciesBuffer[spIndex].append(aLoc.writeOutput());
  				speciesBuffer[spIndex].append("," + aLoc.death + ";\n");
  			}
  			
  			/*<-------HGT Sonia Begin ------> */
  			// For each MultiEpiBac that currently hosts at least one plasmid:
  			if(anAgent instanceof MultiEpiBac )
  			{
  				anEpiBac = (MultiEpiBac) anAgent;
  				if (anEpiBac.plasmidHosted.size() > 0)
  				{
  					// Let 
  					spPlasmid[spIndex]++;
  					for (MultiEpisome epiPlasmid : anEpiBac.plasmidHosted)
  						for (String simPlasmid : aSim.plasmidList)
  							if (epiPlasmid.getName().equals(simPlasmid))
  								spPlasmidTypes[spIndex]
  									  [aSim.plasmidList.indexOf(simPlasmid)]++;
  			}
  			if( anEpiBac.plasmidVector.size() > 0 )
  				spConjugEvents[spIndex] += anEpiBac.plasmidVector.size();
  			}
  		/* <------- HGT Sonia End ------> */
  		}
  		
  		for (Species aSpecies : aSim.speciesList)
  		{
  			spIndex = aSpecies.speciesIndex;
  			// Write the agent_state info for this species to file  
  			speciesBuffer[spIndex].append("</species>\n");
  			bufferStateDeath.write(speciesBuffer[spIndex]);
  			
  			// Collate the agent_Sum info.
  			textBuffer = new StringBuffer();
  			textBuffer.append("<species name=\"");
  			textBuffer.append(aSpecies.speciesName);
  			textBuffer.append("\" header=\"population,mass,growthRate");
  			
  			/*<----HGT Sonia begin---->*/
  			if (Simulator.multiEpi)
  			{
  				textBuffer.append(",plasmidBearing");
  				if (plasmidListSize > 0)
  					for (String plName : aSim.plasmidList)
  						textBuffer.append(",").append(plName);
  			 			textBuffer.append(",conjCount");
  			 }
  			 /*<----HGT Sonia end---->*/
  			
  			textBuffer.append("\" >\n");
  			
  			textBuffer.append(spPop[spIndex]).append(",");
  			textBuffer.append(spMass[spIndex]).append(",");
  			textBuffer.append(spGrowth[spIndex]).append(",");
  			
  			/*<----HGT Sonia begin---->*/
  			if(Simulator.multiEpi)
  			{
  				textBuffer.append(",").append(spPlasmid[spIndex]);
  				for (int c=0; c<spPlasmidTypes[spIndex].length; c++)
  					textBuffer.append(",").append(spPlasmidTypes[spIndex][c]);
  				textBuffer.append(",").append(spConjugEvents[spIndex]);
  			}
  			/*<----HGT Sonia end---->*/
  			
  			textBuffer.append(";\n</species>\n");
  			bufferSumDeath.write(textBuffer);
  		}
	}

	/**
	 * \brief Writes information on the shoving grid to the output string
	 * specified.
	 * 
	 * @param outputString	String that will be output in the result file.
	 * @return	String with information about the grid (resolution and size in
	 * I, J, K directions).
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
	 * \brief Writes information on the environment to the result files if
	 * this is a fluctenv simulation.
	 * 
	 * @return	StringBuffer with information about the environment for
	 * writing to the output file.
	 */
	public StringBuffer writeFluctEnvInformation(StringBuffer outputString)
	{
		outputString.append("<Environment");
		outputString.append(" env=\"").append(FluctEnv.envStatus).append("\"");
		outputString.append("/>\n");
		return outputString;
	}
	
	/**
	 * \brief Checks the grid size, dependent on whether this is a biofilm or
	 * chemostat simulation, and sets suitable dimensions.
	 * 
	 * @param aSimulator	The simulation object used to simulate the
	 * conditions specified in the protocol file.
	 * @param root	The agentGrid markup from the XML file.
	 */
	public void checkGridSize(Simulator aSimulator, XMLParser root) 
	{
		/*
		 * Read in the grid resolution from the XML protocol file.
		 * TODO RC - Why are we reading in the resolution from the xml root,
		 * only to overwrite it with info from the domain?!
		 */
		_res = root.getParamDbl("resolution");
		if ( Simulator.isChemostat )
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
			if ( domain.is3D() ) 
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
	 * \brief Creates the shoving grid, padded of LocatedGroups.
	 * 
	 * Create a vectorized array of spatial groups, build their neighbourhood
	 * and store it. Note: Shoving grid is a padded of LocatedGroups.
	 * 
	 * @param aSimulator	The current simulation object that is recreating
	 * the conditions specified in the protocol file.
	 */
	public void createShovGrid(Simulator aSimulator) 
	{
		_grid = new LocatedGroup[_nTotal];	
		for (int index = 0; index < _nTotal; index++)
			_grid[index] = new LocatedGroup(index, this, aSimulator);
		for ( LocatedGroup lg : _grid )
			lg.init();
		/*
		LogFile.writeLogDebug("Debugging Agentcontainer.createShovGrid");
		for ( LocatedGroup lg : _grid )
			lg.countNullNeighbours();
		*/
	}

	/**
	 * \brief Creates the output species and erosion grids.
	 * 
	 * There is one erosion grid, but one grid for each species in the
	 * simulation.
	 * 
	 * @param aSim	The current simulation object that is recreating the
	 * conditions specified in the protocol file.
	 */
	public void createOutputGrid(Simulator aSim) 
	{
		_erosionGrid = ExtraMath.newDoubleArray(_nI + 2, _nJ + 2, _nK + 2);
		_speciesGrid = new SpatialGrid[aSim.speciesList.size()];
		/*
		 * Create a grid for each species in the simulation
		 */
		for (Species aSpecies : aSim.speciesList)
		{
			_speciesGrid[aSpecies.speciesIndex] =
								domain.createGrid(aSpecies.speciesName, 0.0);
		}
	}

	/* ____________________ EROSION & DETACHMENT __________________________ */

	/**
	 * \brief Perform connected volume filtration (connected to bottom) to
	 * determine which agents should be marked for sloughing.
	 */
	protected void markForSloughing()
	{ 
		/*
		 * cvf is true for connected elements, and false for non-connected.
		 */
		Boolean[] cvf = (new ConnectedVolume(_nI, _nJ, _nK)).computeCvf(_grid);
		int numRemoved = 0;
		Double massRemoved = 0.0;
		/*
		 * Mark as detachable all particles in non-valid map positions.
		 */
		for (int index = 0; index < _nTotal; index++)
			// If it's not connected, remove the agents.
			if ( ! cvf[index] && _grid[index].totalMass > 0 )
			{
				numRemoved += _grid[index].group.size();
				massRemoved += _grid[index].totalMass;
				_grid[index].killAll("detachment");
			}

		LogFile.writeLog("Sloughing " + numRemoved + " ("
				+ ExtraMath.toString(massRemoved, false) + " fg)");
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

		double mass = 0.0;
		double tallyVariable = 0.0;
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
					aLoc.particleMass[iComp] *= 1.0 - ratio;

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
	private Double detFunction(double i, double location)
	{
		if (i < 1.0)
			return ExtraMath.sq(_res - (location % _res));
		else
			return ExtraMath.sq(location % _res);
	}

	/**
	 * \brief Performs the discrete removal of agents over the whole _close
	 * group together, instead of considering each Located Group separately. 
	 * 
	 * Note also that here tallyVariable is cumulative: whatever remains rolls
	 * over to the next time-step. This eliminates time-step issues.
	 * 
	 * @param agentGrid	The agent grid.
	 */
	public void removeOnBorder(AgentContainer agentGrid) 
	{
		/*
		 * Find the border points (erosion and sloughing may have changed the
		 * configuration).
		 */
		_levelset.refreshBorder(true, mySim);
		// List of agents to consider for removal.
		LinkedList<LocatedAgent> detGroup = new LinkedList<LocatedAgent>();
		// For groups on _close list:
		for (LocatedGroup borderElem : _levelset.getBorder())
		{
			/*
			 * Tally up tallyVariable, the approximate amount of mass to
			 * remove, by the ratio variable of each border element.
			 */
			borderElem.erosionRatio = SimTimer.getCurrentTimeStep() /
								_grid[borderElem.gridIndex].erosionTime;
			/*
			 *  i.e. ratio =
			 *  (currentTimeStep * detachmentSpeed * numberFreeNeighbours)/
			 *  	(resolution)
			 */
			borderElem.erosionRatio = Math.min(borderElem.erosionRatio, 1.0);
			tallyVariable += borderElem.totalMass * borderElem.erosionRatio;
			// Add them to detGroup.
			for ( LocatedAgent aLoc : borderElem.group )
				detGroup.add(aLoc);
		} // end of: for (LocatedGroup aBorderElement : _levelset.getBorder())
		
		/*
		 * If the tally is smaller than the smallest agent, no point
		 * calculating detachment priorities.
		 */
		Comparator<Object> comp = new LocatedAgent.totalMassComparator();
		LocatedAgent aLoc = Collections.min(detGroup, comp);
		if ( tallyVariable < aLoc.getTotalMass() )
			return;
		
		// Counter of mass so far removed.
		Double mass = 0.0;
		// Counter of agents removed.
		int nDetach = 0;
		// Calculate detPriority for all agents in the _close list.
		comp = new LocatedAgent.detPriorityComparator();
		for (LocatedGroup borderElem : _levelset.getBorder() )
			calcDetPriority(agentGrid, borderElem, borderElem.erosionRatio);
		// aLoc is the most exposed cell.
		aLoc = Collections.max(detGroup, comp);
		while ( aLoc.getTotalMass() < tallyVariable && 
				! detGroup.isEmpty())
		{
			// Remove agents 1 by 1, in order of decreasing detPriority.
			aLoc = Collections.max(detGroup, comp);
			mass += aLoc.getTotalMass();
			tallyVariable -= aLoc.getTotalMass();
			aLoc.die(false);
			aLoc.death = "detachment";
			nDetach++;
			detGroup.remove(aLoc);
		}
		System.out.println("******************************REMOVE ON BORDER**************************************");
		LogFile.writeLog("Eroding " + nDetach + " ("
				+ ExtraMath.toString(mass, true) + "/"
				+ ExtraMath.toString(tallyVariable, true) + " fg) from "
				+ _levelset.getBorder().size() +" elements.");
	}


	/**
	 * \brief Calculate detachment priorities.
	 * 
	 * @param agentGrid	This agent grid of located agents.
	 * @param aBorderElement	Group of located agents that are on the border.
	 * @param ratio	Ratio value set for that LocatedGroup.
	 */
	public void calcDetPriority(AgentContainer agentGrid,
									LocatedGroup aBorderElement, Double ratio)
	{
		int i = 0;
		// Reset all detPriority values to zero
		for (LocatedAgent aLoc:aBorderElement.group)
			aLoc.detPriority = 0.0;
		/*
		 * For each free neighbour run through the agents in our border
		 * element, adding the square of the agent's proximity to that
		 * neighbour (i.e. distance from opposite neighbour) to its
		 * detPriority variable.
		 */
		for ( i=0; i < 3; i += 2)
		{
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
	 * \brief Check that the coordinates expressed in a ContinuousVector are
	 * defined on this grid.
	 * 
	 * The grid is padded but dc uses shifted coordinates.
	 * 
	 * @return Boolean stating whether or not the ContinuousVector are valid.
	 */
	public Boolean isValid(ContinuousVector cC)
	{
		return isValid(getGridPosition(cC));
	}



	/**
	 * \brief Find the voxel a continuous coordinate lies in and return the
	 * index.
	 * 
	 * @param position	The continuous coordinate to find the index for.
	 * @return Index on the 1D (vectorized) array.
	 */
	public int getIndexedPosition(ContinuousVector position) 
	{
		if ( Simulator.isChemostat )
			return 0;
		else
		{
			//TODO Rob 13Mar2015: Check padding is being dealt with properly.
			int i = (int) Math.floor(position.x / _res) + 1;
			int j = (int) Math.floor(position.y / _res) + 1;
			int k = (int) Math.floor(position.z / _res) + 1;
			return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
			/*
			TODO Rob: started to unify all the ways of converting continuous
			to discrete coordinates but ran into difficulties when values = 0
			
			DiscreteVector dV = new DiscreteVector(position, _res);
			return dV.i + dV.j*(_nI+2) + dV.k*(_nI+2)*(_nJ+2);
			*/
		}
	}

	/**
	 * \brief Find the voxel a discrete vector coordinate lies in and return
	 * the index.
	 * 
	 * @param coord	The discrete vector coordinate to find the index for.
	 * @return Index on the 1D (vectorized) array.
	 */
	public int getIndexedPosition(DiscreteVector coord)
	{
		int i = coord.i + 1;
		int j = coord.j + 1;
		int k = coord.k + 1;

		return i + j * (_nI + 2) + k * (_nI + 2) * (_nJ + 2);
	}

	/**
	 * \brief Takes a ContinuousVector as input (in microns) and returns a
	 * DiscreteVector expressing its location on the grid (in voxels). 
	 * 
	 * TODO Rob 13Mar2015: Compare Spatialgrid.getDiscreteCoordinates()
	 * 
	 * @param position	ContinuousVector that expresses an agent location.
	 * @return	DiscreteVector expressing the agents position in terms of the
	 * agent grid.
	 */
	public DiscreteVector getGridPosition(ContinuousVector position)
	{
		int i = (int) Math.floor(position.x / _res);
		int j = (int) Math.floor(position.y / _res);
		int k = (int) Math.floor(position.z / _res);

		return new DiscreteVector(i, j, k);
	}

	/**
	 * \brief Takes an voxel integer index and returns a DiscreteVector
	 * containing the X,Y, and Z coordinates of that voxel.
	 *  
	 * @param index	Integer index specifying a voxel grid space on the agent
	 * grid.
	 * @return	A discrete vector of the coordinates of this grid.
	 */
	public DiscreteVector getGridPosition(int index)
	{
		//LogFile.writeLogDebug("Debugging AgentContainer.getGridPosition");
		// Remember here that "/" means integer division.
		int i, j, k, num, div;
		num = index;
		div = (_nI + 2) * (_nJ + 2);
		k = num / div;
		num = num % div;
		div = (_nI + 2);
		j = num / div;
		i = num % div;
		// TODO Rob 13Mar2015: Check padding.
		DiscreteVector out = new DiscreteVector(i - 1, j - 1, k - 1);
		//LogFile.writeLogDebug("\tIndex "+index+" -> "+out.toString());
		return out;
	}

	/**
	 * \brief Takes an voxel integer index and returns a ContinuousVector
	 * containing the X, Y, and Z coordinates of the centre of that voxel.
	 * 
	 * @param index	Integer index specifying a voxel grid space on the agent
	 * grid.
	 * @return	A continuous vector of the coordinates of this grid.
	 */
	public ContinuousVector getGridLocation(int index)
	{
		// Remember here that "/" means integer division.
		int i, j, k, num, div;
		num = index;
		div = (_nI + 2) * (_nJ + 2);
		k = num / div;
		num = num % div;
		div = (_nI + 2);
		j = num / div;
		i = num % div;
		// TODO Rob 13Mar2015: Compare ContinuousVector.setToVoxelCenter()
		return new ContinuousVector((i + .5 - 1) * _res, (j + .5 - 1) * _res,
				(k + .5 - 1) * _res);
	}



	/**
	 * \brief Return the agent time step.
	 * 
	 * @return	Double value stating the agent time step.
	 */
	public Double getAgentTimeStep()
	{
		return AGENTTIMESTEP;
	}

	/**
	 * \brief Return the resolution of the grid.
	 * 
	 * @return	Resolution of the agent grid.
	 */
	public Double getResolution()
	{
		return _res;
	}
	
	/**
	 * \brief Return the dimensions of the grid (nI,nJ,nK) as an array.
	 * 
	 * @return	An array containing the dimensions of the grid.
	 */
	public int[] getGridDescription()
	{
		int[] out = { _nI, _nJ, _nK };
		return out;
	}
	
	/**
	 * \brief Return the number of grid cells in the J direction.
	 * 
	 * @return	Number of grid cells in the J direction.
	 */
	public int get_nJ()
	{
		return _nJ;
	}
	
	/**
	 * \brief Return the number of grid cells in the K direction.
	 * 
	 * @return	Number of grid cells in the K direction.
	 */
	public int get_nK()
	{
		return _nK;
	}
	
	/**
	 * \brief Return the shoving grid.
	 * 
	 * @return	LocatedGroup containing the calculated shoving information.
	 */
	
	public LocatedGroup[] getShovingGrid()
	{
		return _grid;
	}
	
	/**
	 * \brief Return the levelset used for modelling detachment.
	 * 
	 * @return	Levelset used for modelling detachment.
	 */
	public LevelSet getLevelSet()
	{
		return _levelset;
	}
	
	/**
	 * \brief Return the status of a given grid voxel, to determine if the
	 * voxel is within a biofilm or not.
	 * 
	 * @param gridVoxel	Index of the grid voxel whose status is being queried.
	 * @return	Integer showing the status of that grid voxel.
	 */
	public int getVoxelStatus(int gridVoxel)
	{
		return _grid[gridVoxel].status;
	}
	
	/**
	 * \brief Return the located group of agents in an agent grid voxel.
	 * 
	 * @param gridVoxel	Integer of the grid voxel for which the located group
	 * is to be returned.
	 * @return	LocatedGroup containing the agents within that grid voxel.
	 */
	public LocatedGroup returnGroupInVoxel(int gridVoxel)
	{
		return _grid[gridVoxel];
	}
}
