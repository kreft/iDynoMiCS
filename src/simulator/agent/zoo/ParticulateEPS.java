/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 */

/**
 * ______________________________________________________
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;

import simulator.agent.LocatedAgent;
import simulator.Simulator;
import simulator.geometry.ContinuousVector;

import utils.ExtraMath;
import utils.XMLParser;

public class ParticulateEPS extends LocatedAgent {

	// Serial version used for the serialisation of the class
	private static final long   serialVersionUID = 1L;

	//private static StringBuffer tempString;

	/* __________________ CONSTRUCTOR _________________________________ */

	/**
	 * Default constructor, only called to create the progenitor of a species
	 */
	public ParticulateEPS() {
		super();
		_speciesParam = new ParticulateEPSParam();
	}

	public Object clone() throws CloneNotSupportedException {
		ParticulateEPS o = (ParticulateEPS) super.clone();
		return o;
	}

	/**
	 * Initialisation procedure based on XML parameter file Used when creating a
	 * progenitor
	 */
	public void initFromProtocolFile(Simulator aSimulator, XMLParser aSpeciesRoot) {
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		_agentGrid = aSimulator.agentGrid;
		init();
	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// this writes no unique values, so doesn't need unique reading-in
		super.initFromResultFile(aSim,singleAgentData);
	}

	public void init() {
		// Lineage management : this is a new agent, he has no known parents
		_generation = 0;
		_genealogy = 0;

		// Determine the radius, volume and total mass of the agent
		updateSize();
	}

	/**
	 * Called by ParticulateEps.createAgent and to obtain another instance of
	 * the same species (totally independent) The returned agent is NOT
	 * registered
	 */
	public ParticulateEPS sendNewAgent() throws CloneNotSupportedException {
		ParticulateEPS baby = (ParticulateEPS) this.clone();
		init();

		return baby;
	}

	/**
	 * Create an agent (who a priori is registered in at least one container;
	 * this agent is located !
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
	 * Called by a capsuled agent when excreting its capsule
	 * @param position
	 * @param mass
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
	
	public boolean willDie() {
		if (_totalMass<0) return true;
		return getRadius(true)<=ExtraMath.deviateFrom(getSpeciesParam().deathRadius,
		        getSpeciesParam().deathRadiusCV);
	}

	public void mutatePop() {
		super.mutatePop();
	}

	public void internalStep() {
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

	public boolean willDivide() {
		return getRadius(true)>getSpeciesParam().divRadius;
	}

	public boolean willTransfer() {
		return getRadius(true)<=ExtraMath.deviateFrom(getSpeciesParam().transferRadius,
		        getSpeciesParam().deathRadiusCV);
	}

	/**
	 * When becoming too small, try to transfer your biomass to your neighbour
	 */
	protected void transferBiomass() {
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
	 * Test if the agent has to die and remove it from any container
	 */
	public void die(boolean isStarving) {
		if (isStarving&_totalMass>0) transferBiomass();		
		super.die(isStarving);
	}

	public void updateVolume() {
		_volume = 0;
		for (int i = 0; i<particleMass.length-1; i++)
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];

		// Add the volume of the EPS capsule to the volume of the intracellular
		// particles
		int i = particleMass.length-1;
		_totalVolume = _volume+particleMass[i]/getSpeciesParam().particleDensity[i];

	}

	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());

		return tempString.toString();
	}

	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());

		return tempString.toString();
	}



	public ParticulateEPSParam getSpeciesParam() {
		return (ParticulateEPSParam) _speciesParam;
	}

	@Override
	protected void conjugate(double elapsedHGTtime) {
		// TODO Auto-generated method stub
		
	}


}
