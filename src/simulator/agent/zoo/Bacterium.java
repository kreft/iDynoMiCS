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

import java.awt.Color;

import org.jdom.Element;

import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

import simulator.Simulator;
import simulator.agent.*;
import simulator.geometry.ContinuousVector;

/**
 * \brief Creates an object of the Bacterium Species. Acts as a superclass for all current species with the exception of ParticulateEPS
 * 
 * Creates an object of the Bacterium Species. Acts as a superclass for all current species with the exception of ParticulateEPS. The 
 * first part of the Bacterium species mark-up in the XML file lists the different particles composing this species; these particles 
 * form ‘compartments’ within the agent that signify the types of biomass that the agent contains. The possible particle types that may 
 * be used are taken from the particle mark-ups defined early in the protocol file, generally consisting of ‘biomass’, ‘inert’, and 
 * ‘capsule’. The initial mass of each particle may be given (in units of femtograms, 1 fg = 10-15 g), but if you set the value to zero 
 * the simulator will assign a reasonable random initial value for the mass. For this Bacterium species, only one particle type is 
 * actually necessary to include, but it is common to include several types. In this species, inert and capsule compartments exhibit 
 * special behaviour: inert biomass will accumulate in the agent, while a capsule will be excreted if it accumulates to take up too 
 * much of the agent’s mass. Note that the ‘capsule’ particle also requires you to specify its type via the class attribute; when the 
 * capsule has become too large and needs to be excreted, an agent of that class will be created (this means that it is necessary to 
 * have previously defined a species with that name.
 * 
 * @author Andreas Dotsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sonia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class Bacterium extends LocatedAgent implements Cloneable 
{	
	/**
	 * Boolean noting whether any EPS particles have been declared as part of this Bacterium
	 */
	protected boolean         _hasEps          = false;
	
	/**
	 * Boolean noting whether any inert particles have been declared as part of this Bacterium
	 */
	protected boolean         _hasInert        = false;
	
	// KA removed 100613 as never referenced
	//protected Species         _inertSpecies;
	
	/**
	 * Previously declared species object of the type of agent that will be created when the capsule is excreted
	 */
	protected Species         _epsSpecies;

	/**
	 * \brief Constructor used to generate progenitor and initialise an object to store relevant parameters
	 * 
	 * Constructor used to generate progenitor and initialise an object to store relevant parameters
	 */
	public Bacterium() 
	{
		super();
		_speciesParam = new BacteriumParam();
	}

	/**
	 * \brief Used by makeKid method to create a daughter Bacterium agent by cloning this agent and parameter objects
	 * 
	 * Used by makeKid method to create a daughter Bacterium agent by cloning this agent and parameter objects
	 * 
	 * @throws CloneNotSupportedException 	Thrown if the agent cannot be cloned
	 */
	public Object clone() throws CloneNotSupportedException {
		Bacterium o = (Bacterium) super.clone();
		o._hasEps = this._hasEps;
		o._epsSpecies = this._epsSpecies;
		return o;
	}

	/**
	 * \brief Creates a Bacterium agent from the parameters specified in the XML protocol file
	 *
	 * Creates a Bacterium agent from the parameters specified in the XML protocol file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) 
	{
		// Initialisation of the Active agent
		super.initFromProtocolFile(aSim, aSpeciesRoot);

		// Check if it is a EPS-producing species
		XMLParser parser;

		for (Element aChild : aSpeciesRoot.getChildren("particle")) {
			// Initialize the xml parser
			parser = new XMLParser(aChild);

			if (parser.getAttributeStr("name").equals("capsule")) {
				_hasEps = true;
				int spIndex = aSim.getSpeciesIndex(parser.getAttributeStr("class"));
				_epsSpecies = aSim.speciesList.get(spIndex);
			}
			if (parser.getAttributeStr("name").equals("inert")) {
				_hasInert = true;

			}						
		}

		/* If no mass defined, use the division radius to find the mass */
		if (this._totalMass==0) {
			guessMass();
			LogFile.writeLog("Guessing "+this.getSpecies().speciesName+" initial mass at: "+this._totalMass);
		}

		// SET CELL RADIUS, GENERATION, AND GENEOLOGY
		init();
	}

	/**
	 * \brief Create an agent using information in a previous state or initialisation file
	 * 
	 * Create an agent using information in a previous state or initialisation file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param singleAgentData	Data from the result or initialisation file that is used to recreate this agent
	 */
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// this writes no unique values, so doesn't need unique reading-in
		// (for a template on how to read in data, look in LocatedAgent.java)
		super.initFromResultFile(aSim,singleAgentData);
	}

	/**
	 * \brief Initialises any new agent (progenitor or daughter cell), setting cell radius, generation, and geneology
	 * 
	 * Initialises any new agent (progenitor or daughter cell), setting cell radius, generation, and geneology
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
	 * \brief Called by Bacterium.createAgent and to obtain another instance of the same species (totally independent). The returned agent is NOT registered
	 * 
	 * Called by Bacterium.createAgent and to obtain another instance of the same species (totally independent). The returned agent is NOT registered
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object cannot be cloned
	 */
	public Bacterium sendNewAgent() throws CloneNotSupportedException 
	{
		// Clone the agent and initialise it
		Bacterium baby = (Bacterium) this.clone();
		baby.init();
		return baby;
	}

	/**
	 * \brief Create a new Bacterium agent (who a priori is registered in at least one container)
	 * 
	 * Create a new Bacterium agent (who a priori is registered in at least one container). This agent is located on the relevant grid
	 */
	public void createNewAgent(ContinuousVector position) 
	{
		try 
		{
			// Get a clone of the progenitor
			Bacterium baby = (Bacterium) sendNewAgent();
			baby.giveName();


			updateMass();
			/* If no mass defined, use the division radius to find the mass */
			// Note this should have been done already in initFromProtocolFile
			if (this._totalMass==0) {
				guessMass();
				LogFile.writeLog("Warning: Bacterium.createNewAgent calling guessMass()");
			}
			
			// randomise its mass
			baby.mutatePop();
			//System.out.println("RADIUS AT THIS POINT: "+this._totalRadius);
			baby.updateSize();
			//System.out.println("RADIUS AFTER CALL: "+this._totalRadius);

			// Just to avoid to be in the carrier
			position.x += this._totalRadius;
			
			baby.setLocation(position);

			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			utils.LogFile.writeLog("Error met in Bacterium:createNewAgent()");
		}
	}

	/**
	 * \brief Mutates inherited parameters and distributes particle mass - either exponentially or normally, dependent on value of distMethod
	 * 
	 *  Mutates inherited parameters and distributes particle mass - either exponentially or normally, dependent on value of distMethod
	 */
	public void mutatePop() 
	{
		// Mutate inherited parameters
		super.mutatePop();

		// distMethod true -> distribute exponentially
		// distMethod false -> distribute normally
		if (getSpeciesParam().distMethod) {
			//LogFile.writeLog("distributing mass exponentially");
			for (int i = 0; i<particleMass.length; i++) {
				particleMass[i] *= ExtraMath.getExp2Rand();
			}
		} else {	
			//LogFile.writeLog("distributing mass normally");
			for (int i = 0; i<particleMass.length; i++) {
				particleMass[i] = ExtraMath.deviateFromCV(1.5*particleMass[i], getSpeciesParam().initialMassCV);
			}
		}
	}

	/**
	 * \brief Mutate inherited agent parameters after agent division.
	 * 
	 *  Mutate inherited agent parameters after agent division. 
	 */
	public void mutateAgent() {
		// Mutate inherited parameters
		super.mutateAgent();

		// Now mutate your parameters
	}

	/**
	 * \brief Used by Bacterium.divide() method to create a daughter cell of this agent
	 * 
	 * Used by Bacterium.divide() method to create a daughter cell of this agent
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object cannot be cloned
	 */
	public void makeKid() throws CloneNotSupportedException {
		super.makeKid();
	}

	/* ___________________ STEP METHODS _______________________________ */

	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls
	 */
	
	/**
	 * \brief Called at each time step of the simulation to compute mass growth and update radius, mass, and volume
	 * 
	 * Called at each time step of the simulation (under the control of the method Step of the class Agent) to compute mass growth 
	 * and update radius, mass, and volume. Also determines whether the agent has reached the size at which it must divide, and 
	 * monitors agent death
	 */
	protected void internalStep() {
		// Compute mass growth over all compartments
		grow();

		updateSize();

		// test if the EPS capsule has to be excreted
		manageEPS();

		// Divide if you have to
		if (willDivide()) divide();

		// Die if you have to
		if (willDie()){
			this.death = "tooSmall";
			die(true);
		}
	}

	/**
	 * \brief Converts the agent into particle EPS and inert on agent death
	 * 
	 * Converts the agent into particle EPS and inert on agent death
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has access to any resources
	 */
	public void die(boolean isStarving) {
		super.die(isStarving);
		if (isStarving) {
			if (_hasEps) excreteEPS(1);
		}
	}

	/**
	 * \brief Excretes EPS particle with 75% of initial mass if the EPS capsule is too thick
	 * 
	 * Excretes EPS particle with 75% of initial mass if the EPS capsule is too thick
	 */
	public void manageEPS() {
		if (!_hasEps) { return; }

		// manage excretion
		if (_volume/_totalVolume<(1-getSpeciesParam().epsMax)) {
			double ratio = ExtraMath.getUniRand(.6, .9);
			excreteEPS(ratio);
		}
	}

	/**
	 * \brief Excretes part of the agents bound EPS as slime EPS
	 * 
	 * Excretes part of the agents bound EPS as slime EPS
	 * 
	 * @param ratio	Ratio of EPS that should be excreted
	 */
	public void excreteEPS(double ratio) {
		int indexEPS = this._agentGrid.mySim.getParticleIndex("capsule");

		// Check the mass to excrete exist
		if (particleMass[indexEPS]*ratio==0) return;

		// Create the particle
		ParticulateEPS eps = (ParticulateEPS) _epsSpecies.getProgenitor();
		// If the particle has been sucessfully created, update your size
		if (eps.createByExcretion(this, ratio)) {
			particleMass[indexEPS] *= (1-ratio);
			updateSize();
		}
	}

	/*
	public void excreteInert(double ratio) {
		int indexInert = this._agentGrid.mySim.getParticleIndex("inert");

		// Check the mass to excrete exist
		// if (particleMass[indexInert]*ratio==0) return;

		// Create the particle
		System.out.println("INERT SPECIES NAME: "+_inertSpecies.getProgenitor().getName());
		ParticulateEPS eps = (ParticulateEPS) _inertSpecies.getProgenitor();
		double totalMass = 0;
		if (eps.createInertByExcretion(this, ratio)) {
			for (int index = 0; index<particleMass.length; index++) {
				totalMass += particleMass[indexInert];
				particleMass[indexInert] = 0;
			}
			eps.particleMass[indexInert] = totalMass;
			updateSize();
		}
	}*/

	/**
	 * \brief Determines if this agent has reached either the radius size limit at which it will die, or a state of zero mass
	 * 
	 * Determines if this agent has reached either the radius size limit at which it will die, or a state of zero mass
	 * 
	 * @return Boolean value noting whether the cell will die (true) or not (false)
	 */
	public boolean willDie() {
		updateRadius();
		if (_totalMass<0) return true;
		// Test cell radius
		if (getRadius(false)<=ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
				getSpeciesParam().deathRadiusCV)) return true;

		// Test inert ratio
		/*
		 * if (_hasInert) { int indexInert =
		 * this._agentGrid.mySim.getParticleIndex("inert"); double ratioInert =
		 * particleMass[indexInert]/(particleMass[indexInert]+particleMass[0]);
		 * if(ratioInert>0.7) return true; }
		 */
		return false;
	}

	/**
	 * \brief Sets the mass of the primary particle of the progenitor to half the mass at-division
	 * 
	 * This sets the mass of the primary particle of the progenitor to half the mass at-division (i.e. the average mass after division, 
	 * regardless of babyMassFrac). The mass-at-division is determined by the particle density and the volume-at-division; volume is 
	 * determined by the divRadius, or the divRadius and the z-resolution if in a 2D simulation.
	 */
	public void guessMass(){
		double divVol;
		// We calculate the mass-at-division
		if (Simulator.isChemostat || _agentGrid.is3D) {
			// In chemostats and 3D the cell is spherical
			divVol = ExtraMath.volumeOfASphere(getSpeciesParam().divRadius);
			//LogFile.writeLog("spherical divVol is "+divVol);
		} else {
			 //In 2D the cell is cylindrical
			divVol = ExtraMath.volumeOfACylinder(getSpeciesParam().divRadius,_species.domain.length_Z);
			//LogFile.writeLog("cylindrical divVol is "+divVol);
		}
		this.particleMass[0] = getSpeciesParam().particleDensity[0]*divVol*0.5;
		updateMass();
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
		// (for a template on how to write the header, look in LocatedAgent.java)
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
		// (for a template on how to write data, look in LocatedAgent.java)
		StringBuffer tempString = new StringBuffer(super.writeOutput());

		return tempString.toString();
	}


	/* _______________ RADIUS, MASS AND VOLUME _____________________ */

	/**
	 * \brief Update the volume of this agent by examining the particle density
	 * 
	 * Update the volume of this agent by examining the particle density
	 */
	public void updateVolume() {
		_totalVolume = 0;
		for (int i = 0; i<particleMass.length; i++) {
			_totalVolume += particleMass[i]/getSpeciesParam().particleDensity[i];
		}

		// Compute volume without EPS capsule
		if (_hasEps) {
			int i = particleMass.length-1;
			_volume = _totalVolume-particleMass[i]/getSpeciesParam().particleDensity[i];
		} else {
			_volume = _totalVolume;
		}
	}

	/**
	 * \brief Return the set of parameters that is associated with the object of this species
	 * 
	 * Return the set of parameters that is associated with the object of this species
	 * 
	 * @return Object of BacteriumParam that stores the parameters associated with this species
	 */
	public BacteriumParam getSpeciesParam() {
		return (BacteriumParam) _speciesParam;
	}

	/**
	 * \brief Determine whether this bacterium contains any eps particles (capsule in the protocol file)
	 * 
	 *  Determine whether this bacterium contains any eps particles (capsule in the protocol file)
	 *  
	 *  @return Boolean noting whether this bacterium object contains eps particles
	 */
	public boolean hasEPS() {
		return _hasEps;
	}

	/**
	 * \brief Determine whether this bacterium contains any inert particles.
	 * 
	 *  Determine whether this bacterium contains any inert particles.
	 *  
	 *  @return Boolean noting whether this bacterium object contains inert particles
	 */
	public boolean hasInert() {
		return _hasInert;
	}

	/**
	 * \brief Return the simulation time at which this agent was created
	 * 
	 * Return the simulation time at which this agent was created
	 * 
	 * @return	Double noting the simulation time at which this agent was created
	 */
	public double getBirthday(){
		return this._birthday;
	}

	/**
	 * \brief Compute the active fraction of the bacterium, ignoring EPS (i.e. only compare active and inert compartments)
	 * 
	 * Compute the active fraction of the bacterium, ignoring EPS (i.e. only compare active and inert compartments)
	 * 
	 * @return Double noting the fraction of the bacterium that is active
	 */
	public double getActiveFrac() {
		if (!hasInert()) return 1.0;

		int indexActive = this._agentGrid.mySim.getParticleIndex("biomass");
		int indexInert = this._agentGrid.mySim.getParticleIndex("inert");

		double val = particleMass[indexActive]/(particleMass[indexActive] + particleMass[indexInert]); 

		if (Double.isNaN(val)) val = 1.0;

		return val;
	}



	/**
	 * \brief Send the colour associated to the species to the defined EPS capsule (if appropriate)
	 * 
	 * Send the colour associated to the species to the defined EPS capsule (if appropriate)
	 * 
	 * @return Color object that this species of Bacterium has been assigned
	 */
	public Color getColorCapsule() {
		if (_epsSpecies==null) return getSpeciesParam().epsColor;
		else return _epsSpecies.color;
	}


	/**
	 * \brief Used for POV-Ray output, defines the colour that this species of Bacterium has been assigned
	 * 
	 * Used for POV-Ray output, defines the colour that this species of Bacterium has been assigned
	 * 
	 * @return Color object that this species of Bacterium has been assigned
	 */
	public Color getColor() {
		return super.getColor();


	}


	//@Override
	//protected void conjugate(double elapsedHGTtime) {
		// TODO Auto-generated method stub

	//}


}
