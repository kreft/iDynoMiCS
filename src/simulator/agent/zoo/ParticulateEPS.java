/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent.zoo;

import simulator.agent.LocatedAgent;
import simulator.Simulator;
import simulator.geometry.ContinuousVector;

import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Creates an object of the Particulate EPS species
 * 
 * Creates an object of the Particulate EPS species. This represents generic extracellular polymers and contains only the 'capsule' 
 * type. This species often includes a hydrolysis reaction, which converts COD bound in EPS into soluble COD. Note that if a bacterium 
 * species includes a capsule of this species of ParticulateEPS, the EPS species MUST be defined first in the protocol file in order 
 * to avoid an error; because of this the example protocol files included in iDynoMiCS list ParticulateEPS as the first defined species.
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class ParticulateEPS extends LocatedAgent {

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long   serialVersionUID = 1L;

	/**
	 * \brief Default constructor, only called to create the progenitor of this species
	 * 
	 * Default constructor, only called to create the progenitor of this species
	 */
	public ParticulateEPS() {
		super();
		_speciesParam = new ParticulateEPSParam();
	}

	/**
	 * \brief Creates a clone of this ParticulateEPS agent
	 * 
	 * Creates a clone of this ParticulateEPS agent
	 * 
	 * @throws CloneNotSupportedException 	Exception thrown if the agent is a type that cannot be cloned
	 */
	public Object clone() throws CloneNotSupportedException 
	{
		ParticulateEPS o = (ParticulateEPS) super.clone();
		return o;
	}

	/**
	 * \brief Initialises the object by reading the relevant species information from the simulation XML protocol file
	 * 
	 *  Initialises the object by reading the relevant species information from the simulation XML protocol file
	 */
	public void initFromProtocolFile(Simulator aSimulator, XMLParser aSpeciesRoot) {
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		_agentGrid = aSimulator.agentGrid;
		init();
	}

	/**
	 * \brief Create a Particulate EPS agent using information in a previous state or initialisation file
	 * 
	 * Create a Particulate EPS agent using information in a previous state or initialisation file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param singleAgentData	Data from the result or initialisation file that is used to recreate this agent
	 */
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		// this writes no unique values, so doesn't need unique reading-in
		super.initFromResultFile(aSim,singleAgentData);
	}

	/**
	 * \brief Initialise an agent object, setting its generation and genealogy and calculating the agent size
	 * 
	 * Initialise an agent object, setting its generation and genealogy and calculating the agent size
	 */
	public void init() 
	{
		// Lineage management : this is a new agent, he has no known parents
		_generation = 0;
		_genealogy = 0;

		// Determine the radius, volume and total mass of the agent
		updateSize();
	}

	/**
	 * \brief Called by ParticulateEps.createAgent to obtain another instance of the same species (totally independent) 
	 * 
	 * Called by ParticulateEps.createAgent to obtain another instance of the same species (totally independent). The returned agent is 
	 * NOT registered
	 * 
	 * @throws CloneNotSupportedException	Thrown if attempting to clone an agent type that does not implement Cloneable
	 * @return	New instance of the ParticulateEPS species
	 */
	public ParticulateEPS sendNewAgent() throws CloneNotSupportedException {
		ParticulateEPS baby = (ParticulateEPS) this.clone();
		init();

		return baby;
	}

	/**
	 * \brief Create a new Particulate EPS agent (who a priori is registered in at least one container)
	 * 
	 * Create a new Particulate EPS agent (who a priori is registered in at least one container)
	 */
	public void createNewAgent(ContinuousVector position) {
		try {
			ParticulateEPS baby = (ParticulateEPS) sendNewAgent();
			baby.mutatePop();
			baby.setLocation(position);
			baby.updateSize();

			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			utils.LogFile.writeLog("Error met in ParticulateEPS:createNewAgent()");
		}
	}

	/**
	 * \brief Creates a new Particulate EPS agent when a capsuled agent excretes that capsule
	 * 
	 * Creates a new Particulate EPS agent when a capsuled agent excretes that capsule. The birth of this agent is registered in the 
	 * relevant simulation grids
	 * 
	 * @param mother	Bacterium object that is excreting an EPS capsule
	 * @param ratio	Ratio of biomass that is being transferred to the capsule
	 */
	public boolean createByExcretion(Bacterium mother, double ratio) {
		try {
			ParticulateEPS baby = (ParticulateEPS) sendNewAgent();
			baby._movement.reset();
			// randomize its mass
			baby.mutatePop();
			baby.updateSize();

			// Give a location to the new agent and register it on the agent
			// grid
			baby.setLocation(mother.getLocation());
			baby.particleMass[baby.particleMass.length-1] = ratio
			        *mother.particleMass[mother.particleMass.length-1];
			baby.updateSize();

			// Compute movement vector
			baby.setDivisionDirection(baby.getInteractDistance(mother)/2);
			//System.out.println(baby._divisionDirection.toString());
			baby._movement.add(baby._divisionDirection);

			// Register the baby in the pathway guilds and the spatial grid
			boolean sucess = !baby.willDie();
			if (sucess) baby.registerBirth();
			return sucess;

		} catch (CloneNotSupportedException e) {
			return false;
		}
	}

	/**
	 * \brief Creates a new inert particle agent when a capsuled agent excretes that capsule
	 * 
	 * Creates a new inert particle agent when a capsuled agent excretes that capsule. The birth of this agent is registered in the 
	 * relevant simulation grids
	 * 
	 * @param mother	Bacterium object that is excreting an EPS capsule
	 * @param ratio	Ratio of biomass that is being transferred to the capsule
	 */
	public boolean createInertByExcretion(Bacterium mother, double ratio) {
		try {
			ParticulateEPS baby = (ParticulateEPS) sendNewAgent();

			// randomize its mass
			baby.mutatePop();
			baby.updateSize();

			// Give a location to the new agent and register it on the agent
			// grid
			baby.setLocation(mother.getLocation());
			baby.particleMass[baby.particleMass.length-2] = ratio
			        *mother.particleMass[mother.particleMass.length-2];
			baby.updateSize();

			// Compute movement vector
			baby.setDivisionDirection(baby.getInteractDistance(mother));
			baby._movement.add(_divisionDirection);

			// Register the baby in the pathway guilds and the spatial grid
			boolean sucess = !baby.willDie();
			if (sucess) baby.registerBirth();
			return sucess;

		} catch (CloneNotSupportedException e) {
			return false;
		}
	}
	
	/**
	 * \brief Determines if this agent has reached either the radius size limit at which it will die, or a state of zero mass
	 * 
	 * Determines if this agent has reached either the radius size limit at which it will die, or a state of zero mass
	 * 
	 * @return Boolean value noting whether the cell will die (true) or not (false)
	 */
	public boolean willDie() {
		if (_totalMass<0) return true;
		return getRadius(true)<=ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
		        getSpeciesParam().deathRadiusCV);
	}

	/**
	 * \brief Mutate any inherited parameters for this particular agent
	 * 
	 * Mutate any inherited parameters for this particular agent. KA June 2013 - not sure this action is implemented
	 */
	public void mutatePop() 
	{
		super.mutatePop();
	}

	/**
	 * \brief Called at each time step of the simulation to compute mass growth and update radius, mass, and volume
	 * 
	 * Called at each time step of the simulation to compute mass growth and update radius, mass, and volume. Also determines whether 
	 * the agent has reached the size at which it must divide, and monitors agent death
	 */
	public void internalStep() 
	{
		// Compute mass growth over all compartments and update radius, mass and
		// volume
		grow();
		updateSize();

		// Divide if you have to
		if (willDivide()) divide();

		// Die if you have to
		if (willTransfer()) transferBiomass();
		if (willDie()) die(true);
	}

	/**
	 * \brief Determines if this agent has reached the radius size limit at which it will divide
	 * 
	 * Determines if this agent has reached the radius size limit at which it will divide
	 * 
	 * @return Boolean value noting whether the cell will divide (true) or not (false)
	 */
	public boolean willDivide() 
	{
		return getRadius(true)>getSpeciesParam().divRadius;
	}

	/**
	 * \brief Determines if the agent can transfer biomass to a neighbour upon cell death
	 * 
	 * Determines if the agent can transfer biomass to a neighbour upon cell death
	 * 
	 * @return	Boolean noting whether a biomass transfer is possible
	 */
	public boolean willTransfer() {
		return getRadius(true)<=ExtraMath.deviateFromCV(getSpeciesParam().transferRadius,
		        getSpeciesParam().deathRadiusCV);
	}

	/**
	 * \brief Method to transfer biomass to a neighbour should the agent become too small
	 * 
	 * Method to transfer biomass to a neighbour should the agent become too small
	 */
	protected void transferBiomass() 
	{
		// Find a neighbour with the same species in your range
		findCloseSiblings(speciesIndex);

		// Remove to big siblings
		int nNb = _myNeighbors.size();
		for (int iNb = 0; iNb<nNb; iNb++) {
			LocatedAgent aLoc = _myNeighbors.removeFirst();
			if (!aLoc.willDivide()) _myNeighbors.add(aLoc);
		}
		if (_myNeighbors.isEmpty()) return;

		// If other particles are around you, transfer your mass
		nNb = _myNeighbors.size();
		double ratio = 0d;
		for (int iNb = 0; iNb<nNb; iNb++) {
			ratio = nNb-iNb;
			ratio = 1/ratio;
			transferCompounds(_myNeighbors.removeFirst(), ratio);
		}
	}

	/**
	 * \brief Tests if the agent has to die and removes it from any container
	 * 
	 * Tests if the agent has to die and removes it from any container
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has access to any resources
	 */
	public void die(boolean isStarving) {
		if (isStarving&_totalMass>0) transferBiomass();		
		super.die(isStarving);
	}

	/**
	 * \brief Update the volume of this agent by examining the particle density
	 * 
	 * Update the volume of this agent by examining the particle density
	 */
	public void updateVolume() 
	{
		_volume = 0;
		for (int i = 0; i<particleMass.length-1; i++)
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];

		// Add the volume of the EPS capsule to the volume of the intracellular
		// particles
		int i = particleMass.length-1;
		_totalVolume = _volume+particleMass[i]/getSpeciesParam().particleDensity[i];

	}

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	public String sendHeader() 
	{
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());

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

		return tempString.toString();
	}


	/**
	 * \brief Return the set of parameters that is associated with the object of this species
	 * 
	 * Return the set of parameters that is associated with the object of this species
	 * 
	 * @return Object of ParticulateParam that stores the parameters associated with this species
	 */
	public ParticulateEPSParam getSpeciesParam() {
		return (ParticulateEPSParam) _speciesParam;
	}

	//@Override
	//protected void conjugate(double elapsedHGTtime) {
		// TODO Auto-generated method stub
		
	//}


}
