/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation and their participation in relevant reactions
 * 
 * Package of utilities that create and manage agents in the simulation and their participation in relevant reactions. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;

import idyno.SimTimer;

import java.util.Iterator;
import java.util.LinkedList;
import java.awt.Color;

import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

import simulator.*;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import simulator.geometry.boundaryConditions.AllBC;

/**
 * \brief Extends ActiveAgent by adding functionality to control agent grid location, agent shoving, agent death and division, and agent movement
 * 
 * Extends ActiveAgent by adding functionality to control agent grid location, agent shoving, agent death and division, and agent movement. 
 * During each global timestep, agent divisions and agent growth lead to many cases where neighbouring agents will overlap. A relaxation 
 * algorithm is used to find iteratively the new overlap-minimising steady state configuration of agent locations at the end of each 
 * timestep. 
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * @author Rob Clegg (rjc096@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public abstract class LocatedAgent extends ActiveAgent implements Cloneable 
{

	/**
	 * Temporary variable storing the distance between two agents
	 */
	protected static ContinuousVector  _diff              = new ContinuousVector();
	
	/**
	 * Temporary store of the new location this cell will move to
	 */
	protected static ContinuousVector  _newLoc            = new ContinuousVector();

	/**
	 * Radius of this agent
	 */
	protected double                   _radius;
	
	/**
	 * Cell radius including any capsules
	 */
	protected double					_totalRadius;
	
	/**
	 * Radius at which this agent will divide
	 */
	protected double				   _myDivRadius;
	
	/**
	 * Radius at which agent death is triggered
	 */
	protected double					_myDeathRadius;
	
	/**
	 * Volume of this agent
	 */
	protected double                   _volume;
	
	/**
	 * Cell volume including any capsules
	 */
	protected double					_totalVolume;
	
	/**
	 * Agent position - continuous coordinates
	 */
	protected ContinuousVector         _location          = new ContinuousVector();
	
	/**
	 * ContinuousVector noting the move that will be applied to the agents position
	 */
	protected ContinuousVector         _movement          = new ContinuousVector();
	
	/**
	 * Direction in which this cell divides
	 */
	protected ContinuousVector         _divisionDirection = new ContinuousVector();
	
	/**
	 * List of neighbouring agents in this agents vicinity
	 */
	protected LinkedList<LocatedAgent> _myNeighbors       = new LinkedList<LocatedAgent>();

	/**
	 * Index of the agent position on the vectorized grid
	 */
	protected int                      _agentGridIndex;
	
	/**
	 * Boolean noting whether this agent is interacting with a surface (true) or not (false)
	 */
	protected boolean                  _isAttached        = false;

	/**
	 * Detachment priority
	 */
	public double                       detPriority = 0;

	/**
	 * Stores the simulation time since the last division check
	 */
	public double                      _timeSinceLastDivisionCheck = Double.MAX_VALUE;

	/**
	 * Distance based probability from a given neighbour (used in HGT). Sonia 8-12-2010
	 */
	public double					_distProb = 0; 								
	
	/**
	 * Cumulative probability as to whether the plasmid will be transferred 
	 */
	public double					_distCumProb = 0; 	


	/**
	 * \brief Constructor used to generate progenitor and initialise an object to store relevant parameters
	 * 
	 * Constructor used to generate progenitor and initialise an object to store relevant parameters
	 */
	public LocatedAgent() {
		super();
		_speciesParam = new LocatedParam();
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Creates a daughter Located Agent by cloning this agent and parameter objects
	 * 
	 * Creates a daughter Located Agent by cloning this agent and parameter objects
	 * 
	 * @throws CloneNotSupportedException 	Thrown if the agent cannot be cloned
	 */
	public Object clone() throws CloneNotSupportedException {
		LocatedAgent o = (LocatedAgent) super.clone();

		o._location = (ContinuousVector) this._location.clone();
		o._movement = (ContinuousVector) this._movement.clone();
		o._divisionDirection = (ContinuousVector) this._divisionDirection
		.clone();
		o._myNeighbors = (LinkedList<LocatedAgent>) this._myNeighbors.clone();

		o._agentGridIndex = this._agentGridIndex;

		return (Object) o;
	}

	/**
	 * \brief Create a new agent in a specified position
	 * 
	 * Create a new agent in a specified position
	 * 
	 * @param position	Vector stating where this agent should be located
	 */
	public void createNewAgent(ContinuousVector position) 
	{
		try 
		{
			// Get a clone of the progenitor
			LocatedAgent baby = (LocatedAgent) sendNewAgent();
			baby.giveName();

			// randomize its mass
			baby.mutatePop();
			
			baby.updateSize();
			
			this._myDivRadius = getDivRadius();
			baby._myDivRadius = getDivRadius();
			baby._myDeathRadius = getDeathRadius();

			// Just to avoid to be in the carrier
			position.x += this._totalRadius;
			
			baby.setLocation(position);

			baby.registerBirth();

		} 
		catch (CloneNotSupportedException e) 
		{
			utils.LogFile.writeLog("Error met in LocAgent:createNewAgent()");
		}
	}

	/**
	 * \brief Registers a created agent into a respective container. Each agent must be referenced by one such container.
	 *  
	 * Registers a created agent into a respective container. Each agent must be referenced by one such container. In this case, the 
	 * species is registered into the agent grid
	 */
	public void registerBirth() {
		// Register on species and reaction guilds
		super.registerBirth();
	}


	/**
	 * \brief Creates an agent of the specified species and notes the grid in which this is assigned
	 *
	 * Creates an agent of the specified species and notes the grid in which this is assigned
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param xmlMarkUp	A species mark-up within the specified protocol file
	 */
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) 
	{	
		super.initFromProtocolFile(aSim, xmlMarkUp);
		
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
		
	}
	
	/**
	 * \brief Create an agent using information in a previous state or initialisation file
	 * 
	 * Create an agent using information in a previous state or initialisation file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param singleAgentData	Data from the result or initialisation file that is used to recreate this agent
	 */
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		// this routine will read data from the end of the singleAgentData array
		// and then pass the remaining values onto the super class

		// Chemostat "if" added by Sonia 27.10.09
		// Rearranged by Rob 10.01.11

		// find the position to start at by using length and number of values read
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;

		if(Simulator.isChemostat){

			// Rob: this is necessary for the case when biofilm agents in one simulation
			// are transferred into a chemostat for the next.
			_location.set(0, 0, 0);

		}else{

			double newAgentX, newAgentY, newAgentZ;
			newAgentX = Double.parseDouble(singleAgentData[iDataStart]);
			newAgentY = Double.parseDouble(singleAgentData[iDataStart+1]);
			newAgentZ = Double.parseDouble(singleAgentData[iDataStart+2]);
			_location.set(newAgentX, newAgentY, newAgentZ);

		}

		// agent size
		_radius      = Double.parseDouble(singleAgentData[iDataStart+3]);
		_totalRadius = Double.parseDouble(singleAgentData[iDataStart+4]);
		
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();

		// now go up the hierarchy with the rest of the data
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];

		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	
	/**
	 * \brief Called at each time step of the simulation to compute cell growth, update size, and monitor cell death and division
	 * 
	 * Called at each time step of the simulation to compute cell growth, update size, and monitor cell death. Also determines whether 
	 * the agent has reached the size at which it must divide
	 */
	protected void internalStep() {
		// Compute mass growth over all compartments
		grow();

		// Apply this mass growth of all compounds on global radius and mass
		updateSize();

		// Divide if you have to
		if (willDivide())
			divide();

		// Die if you have to
		if (willDie())
			die(true);
	}

	/**
	 * \brief Update the radius of the agent from the current mass (and then the volume) of the agent (EPS included)
	 * 
	 * Update the radius of the agent from the current mass (and then the volume) of the agent (EPS included)
	 */
	public void updateSize() 
	{
		// Update the totalMass field (sum of the particles masses)
		updateMass();
		
		if (_totalMass < 0)
			LogFile.writeLog("Warning: negative mass on agent "+_family+", "+_genealogy);

		// Sum of (particles masses / particles density)
		updateVolume();

		// Compute radius according to the volume
		updateRadius();

		//sonia:chemostat
		if(Simulator.isChemostat){
			//don't do the update of attachment/detachment 

		}else{

			// Check if by chance the agent is close enough to a support to be
			// attached

			updateAttachment();
		}
	}

	/**
	 * \brief Captures cell division by making a clone of this agent using the makeKid method
	 * 
	 * Captures cell division by making a clone of this agent using the makeKid method
	 */
	public void divide() {
		try {
			// Create a daughter cell
			makeKid();
		} catch (CloneNotSupportedException e) {
			LogFile.writeLog("Error met in LocatedAgent.divide()");
		}
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where cell division can be triggered
	 * 
	 * Determines whether or not a cell has reached the radius where cell division can be triggered
	 * 
	 * @return	Boolean stating whether cell division should be triggered (true) or not (false)
	 */
	public boolean willDivide() 
	{
		// this ensures that the checks for when to divide don't occur too often;
		// at most they will occur at the rate of AGENTTIMESTEP
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if (_timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep())
			return false;

		// at this point we will actually check whether to divide
		_timeSinceLastDivisionCheck = 0;

		return getRadius(false) > this._myDivRadius;
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where cell death can be triggered
	 * 
	 * Determines whether or not a cell has reached the radius where cell death can be triggered
	 * 
	 * @return	Boolean stating whether cell death should be triggered (true) or not (false)
	 */
	public boolean willDie() {
		if (_totalMass < 0)
			return true;
		return getRadius(false) <= this._myDeathRadius;
	}

	/**
	 * \brief Kills an agent. Called by detachment and starving test
	 * 
	 * Kills an agent. Called by detachment and starving test
	 */
	public void die(boolean isStarving) {
		super.die(isStarving);
	}

	/**
	 * \brief With it determined that cell division will occur, create a new agent from the existing one
	 * 
	 * With it determined that cell division will occur, create a new agent from the existing one
	 * 
	 * @throws CloneNotSupportedException	Thrown if the agent cannot be cloned
	 */
	public void makeKid() throws CloneNotSupportedException {

		// Create the new instance
		LocatedAgent baby = (LocatedAgent) sendNewAgent();
		// Note that mutateAgent() does nothing yet
		baby.mutateAgent();
		
		this._myDivRadius = getDivRadius();
		baby._myDivRadius = getDivRadius();
		baby._myDeathRadius = getDeathRadius();
		
		// Update the lineage
		recordGenealogy(baby);

		// Share mass of all compounds between two daughter cells and compute
		// new size
		divideCompounds(baby, getBabyMassFrac());
		//sonia:chemostat
		if (Simulator.isChemostat){
			// upon division the daughter cells remain with the coordinates of their progenitor

		}else{
			// Compute movement to apply to both cells
			setDivisionDirection(getInteractDistance(baby)/2);

			// move both daughter cells
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}
		// Now register the agent inside the guilds and the agent grid
		baby.registerBirth();
		baby._netVolumeRate = 0;

		
	}

	/**
	 * \brief On agent division, divides the mass between the old and new agent, at a specified fraction
	 * 
	 * On agent division, divides the mass between the old and new agent, at a specified fraction
	 * 
	 * @param baby	The new agent, which is inheriting mass
	 * @param babyMassFrac	The fraction of this agents mass that should be transferred to the new agent
	 */
	public void divideCompounds(LocatedAgent baby, double babyMassFrac) {
		// Choose the division plan and apply position modifications
		for (int i = 0; i<particleMass.length; i++) {
			baby.particleMass[i] *= babyMassFrac;
			this.particleMass[i] *= 1-babyMassFrac;
		}

		// Update radius, mass, volumes and growth rates
		updateSize();
		baby.updateSize();
		
		updateGrowthRates();
		baby.updateGrowthRates();
	}

	/**
	 * \brief On agent division, transfers EPS between the old and new agent, at a specified ratio
	 * 
	 * On agent division, transfers EPS between the old and new agent, at a specified ratio
	 * 
	 * @param baby	The new agent, which is inheriting mass
	 * @param splitRatio	The ratio of the EPS that should be transferred to the new agent
	 */
	public void transferCompounds(LocatedAgent baby, double splitRatio) {
		// Choose the division plan and apply position modifications
		double m;
		for (int i = 0; i<particleMass.length; i++) {
			m = this.particleMass[i]*splitRatio;
			baby.particleMass[i] += m;
			this.particleMass[i] = this.particleMass[i]-m;
		}

		// Update radius, mass and volumes
		updateSize();
		baby.updateSize();
	}

	/**
	 * \brief Mutate any inherited parameters for a population of agents
	 * 
	 * Mutate any inherited parameters for a population of agents. KA June 2013 - not sure this action is implemented
	 */
	public void mutatePop() 
	{
		// Mutate parameters inherited
		super.mutatePop();
		// Now mutate your parameters
	}

	/**
	 * \brief Set the movement vector that states where to put a newly-created particle
	 * 
	 * Set the movement vector that states where to put a newly-created particle
	 * 
	 * @param distance	Distance between the this agent and the new agent
	 */
	public void setDivisionDirection(double distance) {
		double phi, theta;

		phi = 2*Math.PI*ExtraMath.getUniRandDbl();
		theta = 2*Math.PI*ExtraMath.getUniRandDbl();

		_divisionDirection.x = distance*Math.sin(phi)*Math.cos(theta);
		_divisionDirection.y = distance*Math.sin(phi)*Math.sin(theta);
		_divisionDirection.z =(_agentGrid.is3D ? distance*Math.cos(phi):0);
	}

	/* ______________________ SHOVING ___________________________________ */

	/**
	 * \brief Models a mechanical interaction between two located agents. Implemented by extending classes (LocatedAgent)
	 * 
	 * Models a mechanical interaction between two located agents. Implemented by extending classes (LocatedAgent)
	 * 
	 * @param MUTUAL	Whether movement is shared between two agents or applied only to this one
	 * @param shoveOnly	Boolean noting whether this action is shoving (false) or pulling (shrinking biofilm) (true)
	 * @param seq	Whether the move should be applied immediately or wait until the end of the step
	 * @param gain	Double noting change in position
	 * @return	The move to be applied once the shoving or pull calculations have been performed
	 */
	public double interact(boolean MUTUAL, boolean shoveOnly, boolean seq,
			double gain) {
		boolean willShove = false;

		move();

		// rebuild your neighbourhood
		if (shoveOnly)
			getPotentialShovers(getInteractDistance());
		else
			getPotentialShovers(getInteractDistance() + getShoveRadius());

		Iterator<LocatedAgent> iter = _myNeighbors.iterator();
		while (iter.hasNext()) {
			if (shoveOnly)
				willShove |= addPushMovement(iter.next(), MUTUAL, gain);
			else
				willShove |= addSpringMovement(iter.next(), MUTUAL, gain);

		}
		_myNeighbors.clear();

		// Check interaction with surface
		if (_isAttached&!shoveOnly) {

		}

		willShove = isMoving();

		if (seq)
			return move();
		else
			return 0;
	}

	/**
	 * \brief Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector.
	 * 
	 * Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector. 
	 * Both agents are moved of half the overlapping distance in opposite directions.
	 * 
	 * @param aNeighbour	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether shoving is detected (true) or not (false)
	 */
	public boolean addPushMovement(LocatedAgent aNeighbour, boolean isMutual,
			double gain) {
		double d, distance;

		if (aNeighbour == this)
			return false;

		// Build the escape vector and find the distance between you and your
		// neighbourhood
		d = computeDifferenceVector(_location, aNeighbour._location);

		_diff.normalizeVector();

		// Compute effective cell-cell distance
		distance = getShoveRadius()+aNeighbour.getShoveRadius();
		distance += getSpeciesParam().shoveLimit;
		distance = d-distance;

		/* Apply shoving _________________________________________________ */

		// Compute shoving distance for the current agent
		if (distance<=0) {
			if (isMutual) {
				_diff.times(gain*0.5*Math.abs(distance));
				this._movement.add(_diff);
				aNeighbour._movement.subtract(_diff);
			} else {
				_diff.times(Math.abs(gain*distance));
				this._movement.add(_diff);
			}
			return true;
		} else {
			return false;
		}
	}

	/**
	 * \brief Pulling : The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector.
	 * 
	 * The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector. 
	 * 
	 * @param aNeighbor	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether pulling is detected (true) or not (false)
	 */
	public boolean addSpringMovement(LocatedAgent aNeighbor, boolean isMutual,
			double gain) {
		double d, distance, delta;

		if (aNeighbor == this)
			return false;

		// Build the escape vector and find the distance between you and your
		// neighbourhood
		d = computeDifferenceVector(_location, aNeighbor._location);

		_diff.normalizeVector();

		distance = getShoveRadius()+aNeighbor.getShoveRadius();
		distance += getSpeciesParam().shoveLimit;

		delta = d-distance;
		double lMax = _totalRadius;

		if (delta > 0)
			gain *= Math.exp(-delta * 5 / (lMax));
		if (delta > lMax)
			gain = 0;

		/* Apply shoving _________________________________________________ */

		if (isMutual) {
			_diff.times(-0.5*delta*gain);
			this._movement.add(_diff);
			aNeighbor._movement.subtract(_diff);
		} else {
			_diff.times(-delta*gain);
			this._movement.add(_diff);
		}

		return (_movement.norm()>_radius*gain);
	}

	/**
	 * \brief Computes the shortest distance between this agent and another, stored as ContinuousVectors. This may be around the cyclic boundary
	 * 
	 * Computes the distance between this agent and another, stored as ContinuousVectors. This may be around the cyclic boundary
	 * 
	 * @param me	ContinuousVector stating first agent location
	 * @param him	ContinuousVector stating other agent location
	 * @return the shortest movement vector to go from a to b, take into account the cyclic boundary
	 * @see addOverlapMovement
	 * @see addPullMovement works in 2 and 3D
	 */
	public double computeDifferenceVector(ContinuousVector me,
			ContinuousVector him) {
		double gridLength;

		_diff.x = me.x-him.x;
		// check periodicity in X
		gridLength = _species.domain.length_X;
		if (Math.abs(_diff.x) > .5 * gridLength)
			_diff.x -= Math.signum(_diff.x) * gridLength;

		_diff.y = me.y-him.y;
		// check periodicity in Y
		gridLength = _species.domain.length_Y;

		if (Math.abs(_diff.y) > .5 * gridLength)
			_diff.y -= Math.signum(_diff.y) * gridLength;

		if (_agentGrid.is3D) {
			_diff.z = me.z-him.z;
			// check periodicity in Z
			gridLength = _species.domain.length_Z;
			if (Math.abs(_diff.z) > .5 * gridLength)
				_diff.z -= Math.signum(_diff.z) * gridLength;

		} else {
			_diff.z = 0;
		}
		double d = Math.sqrt(_diff.x * _diff.x + _diff.y * _diff.y + _diff.z
				* _diff.z);

		if (d==0) {
			d = 1e-2*_radius;
			_diff.alea(_agentGrid.is3D);
		}

		return d;
	}

	/**
	 * \brief Find neighbouring agents in a range around you
	 * 
	 * Find neighbouring agents in a range around you
	 * 
	 * @param radius	The distance to search around the agent location
	 */
	public void getPotentialShovers(double radius) {
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
	}

	/**
	 * \brief Pick a random neighbour from the _myNeigbors collection
	 * 
	 * Pick a random neighbour from the _myNeigbors collection
	 * 
	 * @return	A randomly picked neighbour (LocatedAgent object) from the list of neighbours
	 */
	public LocatedAgent pickNeighbor() {
		if (_myNeighbors.isEmpty())
			return null;
		else
			return _myNeighbors.get(ExtraMath.getUniRandInt(0, _myNeighbors
					.size()));
	}

	/**
	 * \brief Find a sibling of this agent
	 * 
	 * Find a sibling of this agent
	 * 
	 * @param indexSpecies	The index used to reference this species in the simulation dictionary
	 */
	public void findCloseSiblings(int indexSpecies) 
	{
		int nNb;
		boolean test;
		double shoveDist;
		LocatedAgent aNb;

		getPotentialShovers(getInteractDistance());
		nNb = _myNeighbors.size();

		for (int iNb = 0; iNb<nNb; iNb++) {
			aNb = _myNeighbors.removeFirst();
			// test EPS-species
			test = (indexSpecies==aNb.speciesIndex);

			// Test distance
			shoveDist = 2*(getShoveRadius()+aNb.getShoveRadius());
			test = test
			&& computeDifferenceVector(_location, aNb.getLocation()) <= shoveDist;

			if (test & aNb != this)
				_myNeighbors.addLast(aNb);
		}
	}

	/**
	 * \brief With the agent move calculated, apply this movement, taking care to respect boundary conditions
	 * 
	 * With the agent move calculated, apply this movement, taking care to respect boundary conditions
	 * 
	 */
	public double move() {
		if (!_movement.isValid()) {
			LogFile.writeLog("Incorrect movement coordinates");
			_movement.reset();
		}

		if (!_agentGrid.is3D&&_movement.z!=0) {
			_movement.z = 0;
			_movement.reset();
			LogFile.writeLog("Try to move in z direction !");
		}

		// No movement planned, finish here
		if (_movement.isZero())
			return 0;

		// Test the boundaries
		checkBoundaries();

		// Now apply the movement
		_location.set(_newLoc);
		_agentGrid.registerMove(this);

		double delta = _movement.norm();
		_movement.reset();

		return delta/_totalRadius;
	}

	/**
	 * \brief Used by the move method to determine if an agents move crosses any of the domain's boundaries
	 * 
	 * Used by the move method to determine if an agents move crosses any of the domain's boundaries
	 */
	public void checkBoundaries() 
	{
		// Search a boundary which will be crossed
		_newLoc.set(_location);
		_newLoc.add(_movement);
		AllBC aBoundary = getDomain().testCrossedBoundary(_newLoc);
		int nDim = (_agentGrid.is3D ? 3 : 2);
		boolean boundaryCrossed = (aBoundary!=null);
		int counter = 0;

		// Test all boundaries and apply corrections according to crossed
		// boundaries
		while (boundaryCrossed) {
			counter++;
			aBoundary.applyBoundary(this, _newLoc);
			aBoundary = getDomain().testCrossedBoundary(_newLoc);

            boundaryCrossed = (aBoundary!=null)|(counter>nDim);
			if (counter > nDim)
				System.out.println("LocatedAgent.move() : problem!");
		}
	}

	/**
	 * \brief Mutate inherited agent parameters after agent division.
	 * 
	 * Mutation Function. If you don't want apply a mutation in a specified class, do not redefine this method. If you want, you are 
	 * free to choose which fields to mutate for each different class by a simple redefinition
	 * 
	 */
	public void mutateAgent() {
		// Mutate parameters inherited
		super.mutateAgent();
		// Now mutate your parameters
	}

	/**
	 * \brief Add the reacting concentration of an agent to the received grid
	 * 
	 * Add the reacting concentration of an agent to the received grid
	 * 
	 * @param aSpG	Spatial grid used to sum catalysing mass
	 * @param catalystIndex	Index of the compartment of the cell supporting the reaction
	 */
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) {
		if (isDead)
			return;

		double density = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if (Double.isNaN(density) | Double.isInfinite(density))
			density = 0;
		aSpG.addValueAt(density, _location);
	}

	/**
	 * \brief Add the total concentration of an agent on received grid
	 * 
	 * Add the total concentration of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum total mass
	 */
	public void fitMassOnGrid(SpatialGrid aSpG) 
	{
		if (isDead)
			return;

		double density = _totalMass/aSpG.getVoxelVolume();
		if (Double.isNaN(density) | Double.isInfinite(density))
			density = 0;
		aSpG.addValueAt(density, _location);
	}

	/**
	 * \brief Add the total volume rate of an agent on received grid
	 * 
	 * Add the total volume rate of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum volume
	 */
	public void fitVolRateOnGrid(SpatialGrid aSpG) {
		double volRate = _netVolumeRate/aSpG.getVoxelVolume();
		if (Double.isNaN(volRate) | Double.isInfinite(volRate))
			volRate = 0;
		aSpG.addValueAt(volRate, _location);
	}

	/**
	 * \brief Add the reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * Add the total reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * @param aRateGrid	Spatial grid used to store total reaction rate
	 * @param reactionIndex	Index of this declared reaction in the simulation dictionary
	 */
	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
		if (isDead)
			return;

		// growthRate is in [fgX.hr-1] so convert to concentration:
		// [fgX.um-3.hr-1 = gX.L-1.hr-1]
		double volRate = growthRate[reactionIndex]/aRateGrid.getVoxelVolume();

		if (Double.isNaN(volRate) | Double.isInfinite(volRate))
			volRate = 0;

		aRateGrid.addValueAt(volRate, _location);
	}

	/* _______________ FILE OUTPUT _____________________ */

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",");

		// location info and radius
		tempString.append("locationX,locationY,locationZ,radius,totalRadius");

		return tempString.toString();
	}

	/**
	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * @return	String containing results associated with this agent
	 */
	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());
		tempString.append(",");

		// location info and radius
		tempString.append(_location.x+","+_location.y+","+_location.z+",");
		tempString.append(_radius+","+_totalRadius);

		return tempString.toString();
	}

	/* _______________ RADIUS, MASS AND VOLUME _____________________ */

	/**
	 * \brief Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 * 
	 * Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 */
	public void updateVolume() {
		_volume = 0;
		for (int i = 0; i<particleMass.length; i++) {
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];
		}
		_totalVolume = _volume;
	}

	/**
	 * \brief Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 * 
	 * Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 */
	public void updateRadius() {

		//sonia:chemostat 22.02.2010
		if(Simulator.isChemostat || _species.domain.is3D){
			_radius = ExtraMath.radiusOfASphere(_volume);
			_totalRadius = ExtraMath.radiusOfASphere(_totalVolume);
		}else{
			_radius = ExtraMath.radiusOfACylinder(_volume,
					_species.domain.length_Z);
			_totalRadius = ExtraMath.radiusOfACylinder(_totalVolume,
					_species.domain.length_Z);
		}
	}

	/**
	 * \brief Update the attachment, determining if an agent location crosses any boundaries
	 * 
	 * Update the attachment, determining if an agent location crosses any boundaries	 
	 *  
	 * @return	Boundary that has been crossed
	 */
	public AllBC updateAttachment() {
		// Search a boundary which will be crossed
		double distance;
		for (AllBC aBoundary : getDomain().getAllBoundaries()) {
			if (aBoundary.isSupport()) {
				distance = aBoundary.getDistance(this._location);
				_isAttached = distance<=(3*this._totalRadius);
				return aBoundary;
			}
		}
		return null;
	}

	/**
	 * \brief Add movement to the ContinuousVector storing the agents move
	 * 
	 * Add movement to the ContinuousVector storing the agents move
	 * 
	 * @param aMove	ContinuousVector to add to the movement vector
	 */
	public void addMovement(ContinuousVector aMove) {
		this._movement.add(aMove);
	}

	/**
	 * \brief Return the set of parameters associated with this agent (LocatedParam object)
	 * 
	 * Return the set of parameters associated with this agent (LocatedParam object)
	 * 
	 * @return LocatedParam object of parameters associated with this agent
	 */
	public LocatedParam getSpeciesParam() {
		return (LocatedParam) _speciesParam;
	}

	/**
	 * \brief Return the volume of this agent, with or without the capsule
	 * 
	 *  Return the volume of this agent, with or without the capsule 
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be included in this calculation
	 * @return	Double specifying the volume of this agent
	 */
	public double getVolume(boolean withCapsule) {
		return (withCapsule ? _totalVolume : _volume);
	}

	/**
	 * \brief Return the radius of this agent, with or without the capsule
	 * 
	 *  Return the radius of this agent, with or without the capsule 
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be included in this calculation
	 * @return	Double specifying the radius of this agent
	 */
	public double getRadius(boolean withCapsule) {
		return (withCapsule ? _totalRadius : _radius);
	}

	/**
	 * \brief Return the mass of this agent, with or without the capsule
	 * 
	 *  Return the mass of this agent, with or without the capsule 
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be included in this calculation
	 * @return	Double specifying the mass of this agent
	 */
	public double getMass(boolean withCapsule) {
		return (withCapsule ? _totalMass : _totalMass);
	}

	public double getMaximumRadius() {
		return getSpeciesParam().divRadius
		* (1 + getSpeciesParam().divRadiusCV);
	}

	/**
	 * \brief Determine whether this cell has any EPS
	 * 
	 * Determine whether this cell has any EPS
	 * 
	 * @return	Boolean noting whether this cell has any EPS
	 */
	public boolean hasEPS() 
	{
		return false;
	}

	/**
	 * \brief Determine whether this agent contains any inert particles
	 * 
	 * Determine whether this agent contains any inert particles
	 * 
	 * @return	Boolean noting whether this agent contains any inert particles
	 */
	public boolean hasInert() {
		return false;
	}

	/**
	 * \brief Return the shove factor to be used in shoving for this species of agent
	 * 
	 * Return the shove factor to be used in shoving for this species of agent
	 * 
	 * @return	Double specifying the shove factor that will be applied
	 */
	public double getShoveFactor() {
		return ((LocatedParam) _speciesParam).shoveFactor;
	}

	/**
	 * \brief Return the shove radius to be used in shoving for this species of agent
	 * 
	 * Return the shove radius to be used in shoving for this species of agent
	 * 
	 * @return	Double specifying the shove radius that will be applied
	 */
	public double getShoveRadius() {
		return _totalRadius*((LocatedParam) _speciesParam).shoveFactor;
	}

	/**
	 * \brief Return the shoving interaction distance to be used in shoving for this species of agent
	 * 
	 * Return the shoving interaction distance to be used in shoving for this species of agent
	 * 
	 * @return	Double specifying the shoving interaction distance that will be applied
	 */
	public double getInteractDistance() {
		return 2*getShoveRadius()+((LocatedParam) _speciesParam).shoveLimit;
	}

	/**
	 * \brief Return the shoving interaction distance to be used in shoving against a specified agent
	 * 
	 * Return the shoving interaction distance to be used in shoving against a specified agent
	 * 
	 * @return	Double specifying the shoving interaction distance that will be applied
	 */
	public double getInteractDistance(LocatedAgent baby) {
		return getShoveRadius() + baby.getShoveRadius()
		+ ((LocatedParam) _speciesParam).shoveLimit;
	}

	/**
	 * \brief Return the fraction of mass that is transferred to the new agent on cell division
	 * 
	 * Return the fraction of mass that is transferred to the new agent on cell division
	 * 
	 * @return	Double stating the fraction of mass that is transferred to the new agent on cell division
	 */
	public double getBabyMassFrac() {
		return ExtraMath.deviateFromCV(getSpeciesParam().babyMassFrac,
										getSpeciesParam().babyMassFracCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell division is triggered
	 * 
	 * Return the agent radius at which cell division is triggered
	 * 
	 * @return	Double stating the agent radius at which cell division is triggered
	 */
	public double getDivRadius()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().divRadius,
											getSpeciesParam().divRadiusCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell death is triggered
	 * 
	 * Return the agent radius at which cell death is triggered
	 * 
	 * @return	Double stating the agent radius at which cell death is triggered
	 */
	public double getDeathRadius() {
		return ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
										getSpeciesParam().deathRadiusCV);
	}

	/**
	 * \brief Determine if an agent has a move to perform
	 * 
	 * Determine if an agent has a move to perform
	 * 
	 * @return Boolean noting whether the agent has a move to perform
	 */
	public boolean isMoving() {
		return (_movement.norm()>_totalRadius/10);
	}

	/**
	 * \brief Determine if an agent is attached to a surface
	 * 
	 * Determine if an agent is attached to a surface
	 * 
	 * @return Boolean noting whether the agent is attached to a surface
	 */
	public boolean isAttached() {
		return _isAttached;
	}

	/**
	 * \brief Return the active fraction of this agent
	 * 
	 * Return the active fraction of this agent
	 * 
	 * @return	Double value stating the active fraction of this agent
	 */
	public double getActiveFrac() {
		return 1.0;
	}

	/**
	 * \brief Return the colour assigned to this agent in POV-Ray output
	 * 
	 * Return the colour assigned to this agent in POV-Ray output
	 * 
	 * @return	Colour assigned to this agent as specified in the protocol file
	 */
	public Color getColor() {
		return _species.color;
	}

	/**
	 * \brief Return the colour assigned to any capsules contained in this agent in POV-Ray output
	 * 
	 * Return the colour assigned to any capsules contained in this agent in POV-Ray output
	 * 
	 * @return	Colour assigned to this agent capsules as specified in the protocol file
	 */
	public Color getColorCapsule() {
		return Color.green;
	}

	/**
	 * \brief Return the location of this agent
	 * 
	 * Return the location of this agent
	 * 
	 * @return	ContinuousVector stating the location of this agent
	 */
	public ContinuousVector getLocation() {
		return _location;
	}

	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class detPriorityComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			return (((LocatedAgent) b1).detPriority>((LocatedAgent) b2).detPriority ? 1 : -1);
		}
	}

	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class totalMassComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			return (((LocatedAgent) b1)._totalMass>((LocatedAgent) b2)._totalMass ? 1 : -1);
		}
	}

	/**
	 * \brief Return the distance between two agents
	 * 
	 * Return the distance between two agents
	 * 
	 * @param aLoc	LocatedAgent to find distance to
	 * @return distance between two agents (assuming cyclic boundaries)
	 */
	public double getDistance(LocatedAgent aLoc) {
		return computeDifferenceVector(_location, aLoc._location);
	}

	/**
	 * \brief Set the location of this agent to the supplied continuous vector
	 * 
	 * Set the location of this agent to the supplied continuous vector
	 * 
	 * @param cc	Location which this agent should be assigned to
	 */
	public void setLocation(ContinuousVector cc) 
	{

		//sonia:chemostat
		//set the location of the newborns to zero

		if(Simulator.isChemostat){
			cc.set(0,0,0);
			_location.x = cc.x;
			_location.y = cc.y;
			_location.z = cc.z;

		}else{
			_location.x = cc.x;
			_location.y = cc.y;
			_location.z = cc.z;
		}
	}

	/**
	 * \brief Return the continuous vector that states this agents move
	 * 
	 * Return the continuous vector that states this agents move
	 * 
	 * @return continuous vector that states this agents move
	 */
	public ContinuousVector getMovement() {
		return _movement;
	}

	/**
	 * \brief Return the index of the grid on which this agent is placed
	 * 
	 * Return the index of the grid on which this agent is placed
	 * 
	 * @return Integer grid index of where this agent is placed
	 */
	public int getGridIndex() {
		return _agentGridIndex;
	}

	/**
	 * \brief Return the LocatedGroup of agents that are present in the location where this agent is placed
	 * 
	 * Return the LocatedGroup of agents that are present in the location where this agent is placed
	 * 
	 * @return	LocatedGroup containing all agents present in the same grid space as this agent
	 */
	public LocatedGroup getGridElement() {
		return _agentGrid.getShovingGrid()[_agentGridIndex];
	}

	/**
	 * \brief Move this agent to another grid index
	 * 
	 * Move this agent to another grid index
	 * 
	 * @param aGridIndex	Grid index in which this agent should now be placed
	 */
	public void setGridIndex(int aGridIndex) {
		_agentGridIndex = aGridIndex;
	}

	/**
	 * \brief Return the domain where this agent is contained
	 * 
	 * Return the domain where this agent is contained
	 * 
	 * @return The domain where this agent is contained (Domain object)
	 */
	public Domain getDomain() {
		return _species.domain;
	}

}