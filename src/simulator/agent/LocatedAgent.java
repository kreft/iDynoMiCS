/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

/**
 * ______________________________________________________
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
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

public abstract class LocatedAgent extends ActiveAgent implements Cloneable {

	/* Temporary variables stored in static fields __________________________ */
	protected static ContinuousVector  _diff              = new ContinuousVector();
	protected static ContinuousVector  _newLoc            = new ContinuousVector();

	/* Parameters specific to the agent _____________________________________ */
	protected double                   _radius, _totalRadius;
	protected double				   _myDivRadius, _myDeathRadius;
	protected double                   _volume, _totalVolume;

	/* Agent's location ____________________________________________________ */
	// Agent position and agent movement are expressed with continuous
	// coordinates
	protected ContinuousVector         _location          = new ContinuousVector();
	protected ContinuousVector         _movement          = new ContinuousVector();
	protected ContinuousVector         _divisionDirection = new ContinuousVector();
	protected LinkedList<LocatedAgent> _myNeighbors       = new LinkedList<LocatedAgent>();

	// Index of the agent position on the vectorized grid
	protected int                      _agentGridIndex;
	protected boolean                  _isAttached        = false;

	// Detachment priority
	public double                       detPriority = 0;

	// for timestep issues
	public double                      _timeSinceLastDivisionCheck = Double.MAX_VALUE;

	//sonia 8-12-2010
	// distance based probability from a given neighbour (used in HGT)
	public double					_distProb = 0; 								
	public double					_distCumProb = 0; 	


	/* _______________________ CONSTRUCTOR _________________________________ */

	/**
	 * Empty constructor
	 */
	public LocatedAgent() {
		super();
		_speciesParam = new LocatedParam();
	}

	@SuppressWarnings("unchecked")
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
	 * Create a new agent with mutated parameters based on species default
	 * values and specifies its position
	 */
	/**
	 * Create an agent (who a priori is registered in at least one container;
	 * this agent is located !
	 */
	public void createNewAgent(ContinuousVector position) {
		try {
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

		} catch (CloneNotSupportedException e) {
			utils.LogFile.writeLog("Error met in LocAgent:createNewAgent()");
		}
	}

	/**
	 * Register the agent on the agent grid and on the guilds
	 */
	public void registerBirth() {
		// Register on species and reaction guilds
		super.registerBirth();
	}

	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) {
		super.initFromProtocolFile(aSim, xmlMarkUp);
		
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
		
	}
	
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
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

	/* _____________________HIGH-LEVEL METHODS _____________________________ */

	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls
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
	 * Update the radius of the agent from the current mass (and then the
	 * volume) of the agent (EPS included)
	 */
	public void updateSize() {
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
	 * 
	 */
	public void divide() {
		try {
			// Create a daughter cell
			makeKid();
		} catch (CloneNotSupportedException e) {
			LogFile.writeLog("Error met in LocatedAgent.divide()");
		}
	}

	public boolean willDivide() {
		//jan: commented out since the logic of our simple cell division rule is divide if big enough
		//if (_netGrowthRate<=0) return false;

		// this ensures that the checks for when to divide don't occur too often;
		// at most they will occur at the rate of AGENTTIMESTEP
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if (_timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep())
			return false;

		// at this point we will actually check whether to divide
		_timeSinceLastDivisionCheck = 0;

		return getRadius(false) > this._myDivRadius;
	}

	public boolean willDie() {
		if (_totalMass < 0)
			return true;
		return getRadius(false) <= this._myDeathRadius;
	}

	/**
	 * Kill an agent. Called by detachment and starving test
	 */
	public void die(boolean isStarving) {
		super.die(isStarving);
	}

	/* ________________________________________________________________ */
	/**
	 * Create a new agent from an existing one
	 * 
	 * @throws CloneNotSupportedException
	 *             Called by LocatedAGent.divide()
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

	public void mutatePop() {
		// Mutate parameters inherited
		super.mutatePop();
		// Now mutate your parameters
	}

	/**
	 * Set movement vector to put a new-created particle
	 * 
	 * @param myBaby
	 * @param distance
	 */
	public void setDivisionDirection(double distance) {
		double phi, theta;

		phi = 2*Math.PI*ExtraMath.getUniRand();
		theta = 2*Math.PI*ExtraMath.getUniRand();

		_divisionDirection.x = distance*Math.sin(phi)*Math.cos(theta);
		_divisionDirection.y = distance*Math.sin(phi)*Math.sin(theta);
		_divisionDirection.z =(_agentGrid.is3D ? distance*Math.cos(phi):0);
	}

	/* ______________________ SHOVING ___________________________________ */

	/**
	 * Mechanical interaction between two located agents
	 * 
	 * @param aGroup :
	 *            neighbourhood of the agent
	 * @param MUTUAL :
	 *            movement shared between 2 agents or applied only to this one
	 * @pull : false for shoving, true for pulling (shrinking biofilm)
	 * @seq : apply immediately the movement or waits the end of the step
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
	 * Mutual shoving : The movement by shoving of an agent is calculated based
	 * on the cell overlap and added to the agents movement vector. Both agents
	 * are moved of half the overlapping distance in opposite directions.
	 * 
	 * @param aNeighbour
	 *            reference to the potentially shoving neighbour
	 * @return true, if a shoving is detected
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
	 * 
	 * @param aNeighbor
	 * @param isMutual
	 * @return
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
	 * 
	 * @return
	 */
	public boolean addSpringAttachment() {
		AllBC mySupport = updateAttachment();
		double d, distance, delta;

		d = computeDifferenceVector(_location, mySupport
				.getOrthoProj(_location));
		_diff.normalizeVector();

		distance = _totalRadius*getShoveFactor();
		delta = d-distance;

		/* Apply elastic interaction _______________________________________ */
		double gain = 0.1;
		if (delta < 0)
			gain = 0.1;
		if (delta > 0)
			gain = 0.1;
		if (delta > _totalRadius)
			gain = 0;

		_diff.times(-delta*gain);
		this._movement.add(_diff);

		if (_movement.norm()>_radius*0.1) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @param ContinuousVector
	 *            a location
	 * @param ContinuousVector
	 *            a location
	 * @return the shortest movement vector to go from a to b, take into account
	 * the cyclic boundary
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
	 * Look for neighbours in a range around you
	 */
	public void getPotentialShovers(double radius) {
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
	}

	/**
	 * Pick randomly a Neighbor from the _myNeigbors collection
	 * 
	 * @return
	 */
	public LocatedAgent pickNeighbor() {
		if (_myNeighbors.isEmpty())
			return null;
		else
			return _myNeighbors.get(ExtraMath.getUniRandInt(0, _myNeighbors
					.size()));
	}

	/**
	 * Find a sibling
	 * 
	 * @param indexSpecies
	 * @return
	 */
	public void findCloseSiblings(int indexSpecies) {
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
	 * Apply the movement stored taking care to respect boundary conditions
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

	public void checkBoundaries() {
		// Search a boundary which will be crossed
		_newLoc.set(_location);
		_newLoc.add(_movement);
		AllBC aBoundary = getDomain().testCrossedBoundary(_newLoc);
		int nDim = (_agentGrid.is3D ? 3 : 2);
		boolean test = (aBoundary!=null);
		int counter = 0;

		// Test all boundaries and apply corrections according to crossed
		// boundaries
		while (test) {
			counter++;
			aBoundary.applyBoundary(this, _newLoc);
			aBoundary = getDomain().testCrossedBoundary(_newLoc);

			test = (aBoundary!=null)|(counter>nDim);
			if (counter > nDim)
				System.out.println("LocatedAgent.move() : problem!");
		}
	}

	/* ____________________CELL DIVISION __________________________________ */

	/**
	 * Mutation Function If you don't want apply a mutation in a specified
	 * class, do not redefine this method. If you want, you are free to choose
	 * which fields to mutate for each different class by a simple redefinition
	 * 
	 * @param alea
	 */
	public void mutateAgent() {
		// Mutate parameters inherited
		super.mutateAgent();
		// Now mutate your parameters
	}

	/**
	 * Add the reacting CONCENTRATION of an agent on the received grid
	 * 
	 * @param aSpG :
	 *            grid used to sum catalysing mass
	 * @param catalyst
	 *            index : index of the compartment of the cell supporting the
	 *            reaction
	 */
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) {
		if (isDead)
			return;

		double value = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * Add the total CONCENTRATION of an agent on received grid
	 * 
	 * @param aSpG :
	 *            grid used to sum catalysing mass
	 */
	public void fitMassOnGrid(SpatialGrid aSpG) {
		if (isDead)
			return;

		double value = _totalMass/aSpG.getVoxelVolume();
		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0;
		aSpG.addValueAt(value, _location);
	}

	public void fitVolRateOnGrid(SpatialGrid aSpG) {
		double value;
		value = _netVolumeRate/aSpG.getVoxelVolume();
		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0;
		aSpG.addValueAt(value, _location);
	}

	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
		if (isDead)
			return;

		// growthRate is in [fgX.hr-1] so convert to concentration:
		// [fgX.um-3.hr-1 = gX.L-1.hr-1]
		double value = growthRate[reactionIndex]/aRateGrid.getVoxelVolume();

		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0;

		aRateGrid.addValueAt(value, _location);
	}

	/* _______________ FILE OUTPUT _____________________ */


	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",");

		// location info and radius
		tempString.append("locationX,locationY,locationZ,radius,totalRadius");

		return tempString.toString();
	}

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
	 * Compute the volume on the basis of the mass and density of different
	 * compounds defined in the cell
	 */
	public void updateVolume() {
		_volume = 0;
		for (int i = 0; i<particleMass.length; i++) {
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];
		}
		_totalVolume = _volume;
	}

	/**
	 * Compute the radius on the basis of the volume The radius evolution is
	 * stored in deltaRadius (used for shrinking)
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

	public void addMovement(ContinuousVector aMove) {
		this._movement.add(aMove);
	}

	/* __________________ ACCESSORS ___________________________________ */
	public LocatedParam getSpeciesParam() {
		return (LocatedParam) _speciesParam;
	}

	public double getVolume(boolean withCapsule) {
		return (withCapsule ? _totalVolume : _volume);
	}

	public double getRadius(boolean withCapsule) {
		return (withCapsule ? _totalRadius : _radius);
	}

	public double getMass(boolean withCapsule) {
		return (withCapsule ? _totalMass : _totalMass);
	}

	public double getMaximumRadius() {
		return getSpeciesParam().divRadius
		* (1 + getSpeciesParam().divRadiusCV);
	}

	public boolean hasEPS() {
		return false;
	}

	public boolean hasInert() {
		return false;
	}

	public double getShoveFactor() {
		return ((LocatedParam) _speciesParam).shoveFactor;
	}

	public double getShoveRadius() {
		return _totalRadius*((LocatedParam) _speciesParam).shoveFactor;
	}

	public double getInteractDistance() {
		return 2*getShoveRadius()+((LocatedParam) _speciesParam).shoveLimit;
	}

	public double getInteractDistance(LocatedAgent baby) {
		return getShoveRadius() + baby.getShoveRadius()
		+ ((LocatedParam) _speciesParam).shoveLimit;
	}

	public double getBabyMassFrac() {
		return ExtraMath.deviateFrom(getSpeciesParam().babyMassFrac,
				getSpeciesParam().babyMassFracCV);
	}
	
	public double getDivRadius() {
		return ExtraMath.deviateFrom(getSpeciesParam().divRadius, getSpeciesParam().divRadiusCV);
	}
	
	public double getDeathRadius() {
		return ExtraMath.deviateFrom(getSpeciesParam().deathRadius, getSpeciesParam().deathRadiusCV);
	}

	public boolean isMoving() {
		return (_movement.norm()>_totalRadius/10);
	}

	public boolean isAttached() {
		return _isAttached;
	}

	public double getActiveFrac() {
		return 1.0;
	}

	public Color getColor() {
		return _species.color;
	}

	public Color getColorCapsule() {
		return Color.green;
	}

	public ContinuousVector getLocation() {
		return _location;
	}

	/**
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class detPriorityComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			return (((LocatedAgent) b1).detPriority>((LocatedAgent) b2).detPriority ? 1 : -1);
		}
	}

	/**
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class totalMassComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			return (((LocatedAgent) b1)._totalMass>((LocatedAgent) b2)._totalMass ? 1 : -1);
		}
	}

	/**
	 * @param aLoc
	 * @return distance bw 2 agents assuming cyclic boundaries
	 */
	public double getDistance(LocatedAgent aLoc) {
		return computeDifferenceVector(_location, aLoc._location);
	}

	public void setLocation(ContinuousVector cc) {

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

	public ContinuousVector getMovement() {
		return _movement;
	}

	public int getGridIndex() {
		return _agentGridIndex;
	}

	public LocatedGroup getGridElement() {
		return _agentGrid.getShovingGrid()[_agentGridIndex];
	}

	public void setGridIndex(int aGridIndex) {
		_agentGridIndex = aGridIndex;
	}

	public Domain getDomain() {
		return _species.domain;
	}

}
