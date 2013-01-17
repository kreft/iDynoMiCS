/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */

package simulator.agent;

import idyno.SimTimer;
import java.util.ArrayList;
import org.jdom.Element;

import utils.LogFile;
import utils.XMLParser;
import utils.ExtraMath;
import simulator.Simulator;
import simulator.reaction.Reaction;
import simulator.SpatialGrid;

public abstract class ActiveAgent extends SpecialisedAgent implements HasReaction {

	// Parameters common to all agents of this class

	// Parameters common (strict egality) to all agents of a Species

	// massic growth rate of the agent (the sum of the growth rate of all of
	// its compounds)

	public Reaction[]            allReactions;
	protected ArrayList<Integer> reactionActive;
	protected ArrayList<Integer> reactionKnown;

	protected double             _netGrowthRate = 0;
	protected double             _netVolumeRate = 0;

	protected double[]           growthRate;
	// Reaction parameters : (potentially)mutated from species parameters
	public double[][]            soluteYield;
	public double[][]            reactionKinetic;
	public double[][]            particleYield;

	// Mass of the agent (table for all particles belonging to the agent)
	public double[]              particleMass;
	// Sum of masses of all particles
	protected double             _totalMass;

	/* ________________________ CONSTRUCTOR _________________________________ */

	/**
	 * The constructor is used to create the progenitor of the species
	 */
	public ActiveAgent() {
		super();
		_speciesParam = new ActiveParam();

	}

	/**
	 * Initialize Reaction fields Called by CreateSpecies to initalize the
	 * progenitor The species parameter have already been defined
	 */
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) {
		// Initialisation common to all specialised agents
		super.initFromProtocolFile(aSim, xmlMarkUp);

		/* Create internal compounds________________________________________ */

		// Initialize tables for the compartments description
		int nParticle = aSim.particleDic.size();
		int nReaction = aSim.reactionList.length;
		int nSolute = aSim.soluteList.length;
		int reacIndex;

		particleMass = new double[nParticle];

		// Build the list of particles
		XMLParser parser;
		int particleIndex;

		for (Element aChild : xmlMarkUp.getChildren("particle")) {
			// Initialize the xml parser
			parser = new XMLParser(aChild);
			particleIndex = aSim.getParticleIndex(parser.getAttribute("name"));

			// Set the average mass of the particle within the initial
			// population
			particleMass[particleIndex] = parser.getParamMass("mass");
		}

		updateMass();

		/* Create description of reactions _________________________________ */

		// Initialize the arrays
		allReactions = aSim.reactionList;
		reactionKnown = new ArrayList<Integer>();
		reactionActive = new ArrayList<Integer>();
		growthRate = new double[nReaction];

		soluteYield = new double[nReaction][nSolute];
		reactionKinetic = new double[nReaction][];
		particleYield = new double[nReaction][nParticle];

		// Read the XML file

		for (Element aReactionMarkUp : xmlMarkUp.buildSetMarkUp("reaction")) {
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			Reaction aReaction = allReactions[reacIndex];

			// Add the reaction to the list of known (and active) reactions
			reactionKnown.add(reacIndex);
			if (aReactionMarkUp.getAttributeValue("status").equals("active")) {
				reactionActive.add(reacIndex);
			}

			// If reaction parameters have been redefined, load them ; else load
			// the parameters defined for the reaction
			if (aReactionMarkUp.getContentSize()==0) {
				soluteYield[reacIndex] = aReaction.getSoluteYield();
				particleYield[reacIndex] = aReaction.getParticulateYield();
				reactionKinetic[reacIndex] = aReaction.getKinetic();
			} else {
				aReaction.initFromAgent(this, aSim, new XMLParser(aReactionMarkUp));
			}
		}

		// Now copy these value in the speciesParam strucure
		getSpeciesParam().soluteYield = soluteYield.clone();
		getSpeciesParam().particleYield = particleYield.clone();
		getSpeciesParam().reactionKinetic = reactionKinetic.clone();

	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// this routine will read data from the end of the singleAgentData array
		// and then pass the remaining values onto the super class

		// find the position to start at by using length and number of values read
		int nValsRead = 2 + particleMass.length;
		int iDataStart = singleAgentData.length - nValsRead;

		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT

		// Particle Masses
		for (int iComp = 0; iComp<particleMass.length; iComp++)
			particleMass[iComp] = Double.parseDouble(singleAgentData[iDataStart+iComp]);

		// other values
		_netGrowthRate = Double.parseDouble(singleAgentData[iDataStart+particleMass.length]);
		_netVolumeRate = Double.parseDouble(singleAgentData[iDataStart+particleMass.length+1]);

		// now go up the hierarchy with the rest of the data
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);

		// finally some creation-time calls
		updateSize();
		registerBirth();		
	}	

	public void mutatePop() {
		// Mutate parameters inherited
		super.mutatePop();
		// Now mutate your own class parameters
	}

	/**
	 * Create an agent (who a priori is registered in at least one container;
	 * this agent is NOT located ! Implemented here for compatibility reasons
	 */
	public void createNewAgent() {
		try {
			ActiveAgent baby = (ActiveAgent) sendNewAgent();
			baby.mutatePop();

			// Register the baby in the pathway guilds an
			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			System.out.println("At ActiveAgent: createNewAgent error " + e);
		}
	}

	public void registerBirth() {
		super.registerBirth();
		// register the agent in the metabolic containers
		registerOnAllActiveReaction();
	}

	/* ___________________ DIVISION ______________________________ */

	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException {
		ActiveAgent o = (ActiveAgent) super.clone();
		// Shallow copy : the reaction are not cloned
		o.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
		o.reactionKnown = (ArrayList<Integer>) this.reactionKnown.clone();
		o.allReactions = this.allReactions.clone();

		o.growthRate = new double[growthRate.length];

		o.soluteYield = new double[soluteYield.length][];
		for (int iter = 0; iter<soluteYield.length; iter++) {
			o.soluteYield[iter] = this.soluteYield[iter].clone();
		}

		o.reactionKinetic = new double[reactionKinetic.length][];
		o.particleYield = new double[particleYield.length][];

		for (int iter = 0; iter<reactionKnown.size(); iter++) {
			int jReac = reactionKnown.get(iter);
			if (this.reactionKinetic[jReac]!=null) o.reactionKinetic[jReac] = this.reactionKinetic[jReac]
			                                                                                       .clone();
			o.particleYield[jReac] = this.particleYield[jReac].clone();
		}

		o.particleMass = this.particleMass.clone();

		return (Object) o;
	}

	public void mutateAgent() {
		// Mutate parameters inherited
		super.mutateAgent();
		// Now mutate your own class parameters

	}

	public void die(boolean isStarving) {
		super.die(isStarving);
		// If you are too small, you must die !
		// Decrease the population of your species


		// Unregister from the metabolic guilds
		unregisterFromAllActiveReactions();
	}

	/* ___________________ STEP ______________________________ */
	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls
	 */
	protected void internalStep() {
		grow();
		updateSize();
	}

	/**
	 * Perform growth by calling all active pathways.
	 */
	protected void grow() {

		double deltaMass = 0;
		double[] deltaParticle = new double[particleMass.length];
		int reacIndex = 0;
		double tStep = SimTimer.getCurrentTimeStep();
		double catMass = 0; // Catalyst mass
		double catYield =0;
		_netGrowthRate = 0;
		_netVolumeRate = 0;

		// Compute mass growth rate of each active reaction
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			// Compute the growth rate
			reacIndex = reactionActive.get(iReac);
			catMass = particleMass[allReactions[reacIndex]._catalystIndex];
			// get the growth rate in [fgX.hr-1]
			growthRate[reacIndex] = allReactions[reacIndex].computeSpecGrowthRate(this);

			for (int i = 0; i<particleYield[reacIndex].length; i++) {
				deltaMass = catMass * particleYield[reacIndex][i]*growthRate[reacIndex];
				_netGrowthRate += deltaMass;
				_netVolumeRate += deltaMass/getSpeciesParam().particleDensity[i];
				
				if (allReactions[reacIndex].autocatalytic){
					// Exponential growth/decay
					catYield = particleYield[reacIndex][allReactions[reacIndex]._catalystIndex];
					
					deltaParticle[i] += catMass * (particleYield[reacIndex][i]/catYield)
									    	* Math.expm1(catYield * growthRate[reacIndex]*tStep);
				} else {
					// Constant growth/decay
					deltaParticle[i] += deltaMass*tStep;
				}
			}
		}
		
		// We adjust the particle masses after calculating all the deltaParticle values
		// so that the reactions occur simultaneously
		for (int i = 0; i<particleMass.length; i++){
			particleMass[i] += deltaParticle[i];
		}
		
	}
	
	/**
	 * Take into account all growth
	 */
	public void updateSize() {
		updateMass();
	}

	public void updateMass() {
		_totalMass = ExtraMath.sumVector(particleMass);
	}

	/* ______________________ REACTION MANAGEMENT __________________________ */

	public void addReaction(Reaction aReaction, boolean useDefaultParam) {
		// Add the reaction to the list of known reaction
		reactionKnown.add(aReaction.reactionIndex);

		// Test if specific parameters exist for this reaction
		int index = aReaction.reactionIndex;
		boolean test = getSpeciesParam().soluteYield[index]==null;

		if (useDefaultParam||test) {
			// Use parameters defined in the reaction object
			reactionKinetic[index] = aReaction.getKinetic();
			soluteYield[index] = aReaction.getSoluteYield();
			particleYield[index] = aReaction.getParticulateYield();
		} else {

			// Use parameters defined in the speciesParam structure
			reactionKinetic[index] = getSpeciesParam().reactionKinetic[index];
			soluteYield[index] = getSpeciesParam().soluteYield[index];
			particleYield[index] = getSpeciesParam().particleYield[index];
		}
	}

	public void addActiveReaction(Reaction aReaction, boolean useDefaultParam) {
		addReaction(aReaction, useDefaultParam);
		switchOnReaction(aReaction);
	}

	public void removeReaction(Reaction aPathway) {
		switchOffreaction(aPathway);
		reactionKnown.remove(aPathway);
	}

	// bvm 27.11.08: added the '.reactionIndex' calls to the two
	// lines below in order to get this function to work correctly
	public void switchOffreaction(Reaction aPathway) {
		if (reactionActive.contains(aPathway.reactionIndex)) {
			// need to remove using indexOf because the remove(object) version thinks
			// the int being passed in is the index to remove rather than the object to remove
			reactionActive.remove(reactionActive.indexOf(aPathway.reactionIndex));
			aPathway.removeAgent(this);
		}
	}

	// bvm 27.11.08: added the if statement to prevent adding a reaction if
	// it is already present
	public void switchOnReaction(Reaction aReaction) {
		//		System.out.println("Turn it on? "+aReaction.reactionName);
		if (!reactionActive.contains(aReaction.reactionIndex)) {
			//			System.out.println("Turn on: "+aReaction.reactionName);
			reactionActive.add(aReaction.reactionIndex);
			aReaction.addAgent(this);
		}
	}

	/**
	 * Register the agent on each guild of its activated pathways. Called by
	 * makeKid
	 */
	public void registerOnAllActiveReaction() {
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			allReactions[reactionActive.get(iReac)].addAgent(this);
		}
	}

	/**
	 * Called by the die method
	 */
	public void unregisterFromAllActiveReactions() {
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			allReactions[reactionActive.get(iReac)].removeAgent(this);
		}
	}

	/**
	 * Add the mass of an agent on received grid
	 * @param aSpG : grid used to sum catalysing mass
	 * @param catalyst index : index of the compartment of the cell supporting
	 * the reaction TODO
	 */
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex) {
	}

	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex) {
	}

	/**
	 * Add the mass of an agent on received grid
	 * @param aSpG : grid used to sum catalysing mass TODO
	 */
	public void fitMassOnGrid(SpatialGrid aSpG) {
	}

	/* _______________ FILE OUTPUT _____________________ */

	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",");

		// particle types
		for (int i = 0; i<particleMass.length; i++) {
			tempString.append(_species.currentSimulator.particleDic.get(i));
			tempString.append(",");
		}
		tempString.append("growthRate,volumeRate");
		return tempString.toString();
	}

	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());
		tempString.append(",");

		// Mass of different particles
		for (int i = 0; i<particleMass.length; i++) {
			tempString.append(particleMass[i]);
			tempString.append(",");
		}
		// Agent growth and volume rates
		tempString.append(_netGrowthRate+","+_netVolumeRate);

		return tempString.toString();
	}

	/* ____________________ ACCESSORS & MUTATORS _________________________ */
	public double getTotalMass() {
		return _totalMass;
	}

	public double getParticleMass(int particleIndex) {
		return particleMass[particleIndex];
	}

	public double getNetGrowth() {
		return _netGrowthRate;
	}

	public double getVolGrowth() {
		return _netVolumeRate;
	}

	public void setNetGrowth(double value) {
		_netGrowthRate = value;
	}

	public double[] getSoluteYield(int indexReaction) {
		return soluteYield[indexReaction];
	}

	public double[] getReactionKinetic(int indexReaction) {
		return reactionKinetic[indexReaction];
	}

	public ActiveParam getSpeciesParam() {
		return (ActiveParam) _speciesParam;
	}



}
