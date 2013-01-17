/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.agent;

import idyno.SimTimer;
import simulator.Simulator;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

public abstract class Agent implements Cloneable {

	/* Parameters common to all agents of this class ________________________ */
	/* Temporary variables stored in static fields __________________________ */
	/* Parameters common (strict equality) to all agents of a Species _________ */
	/* Parameters mutated from species parameters ___________________________ */

	/* Parameters specific to the agent _____________________________________ */
	// When has this agent been stepped for the last time ?
	protected int        _lastStep;

	/**
	 * Lineage management
	 * @field _generation is the number of generations between the progenitor
	 * and the current agent,
	 * @field _genealogy is the integer for the binar reading of the 0 and 1
	 * coding the lineage. When a cells divides, one daughter has the index
	 * value 1, the other the index value 0, then this index is added on the
	 * left of the lineage description
	 */
	protected int        _generation = 0;
	protected long        _genealogy  = 0;
	protected int        _family     = 0;
	protected static int nextFamily  = 0;

	protected double     _birthday;

	/* ______________ BUILDER ___________________________________________ */
	// Empty builder, used by example to create a progenitor
	public Agent() {
		_birthday = SimTimer.getCurrentTime();
		_lastStep = SimTimer.getCurrentIter()-1;
	}

	public void initFromProtocolFile(Simulator aSimulator, XMLParser aSpeciesRoot) {

	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT
		_family     = Integer.parseInt(singleAgentData[0]);
		_genealogy  = Long.parseLong(singleAgentData[1]);
		_generation = Integer.parseInt(singleAgentData[2]);
		_birthday   = Double.parseDouble(singleAgentData[3]);

		// (this is top of the hierarchy, so no call to super)
	}

	public void mutateAgent() {
		// Now mutate your parameters
		// TODO
	}

	/**
	 * Create a new agent from an existing one
	 * @throws CloneNotSupportedException
	 */
	public void makeKid() throws CloneNotSupportedException {

		Agent anAgent = (Agent) this.clone();		
		anAgent.mutateAgent();

		// Now register the agent in the appropriate container
		registerBirth();
	}

	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	/**
	 * A created agent has to be referenced by at least one container
	 */
	public abstract void registerBirth();

	/* __________________________ METHODS _______________________________ */
	public void step() {
		_lastStep = SimTimer.getCurrentIter();
		internalStep();
	}

	//sonia 01/2011
	// this was created here so that we can call it during agent step() 
	//in AgentContainer class. The same applies for conjugate()
	public void HGTstep(double elapsedHGTtime){
		conjugate(elapsedHGTtime);
	}

	protected abstract void conjugate(double elapsedHGTtime);
	
	protected abstract void internalStep();



	/* _______________ FILE OUTPUT _____________________ */

	public String sendHeader() {
		// return the header file for this agent's values	
		StringBuffer tempString = new StringBuffer("family,genealogy,generation,birthday");
		return tempString.toString();
	}

	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer("");
		tempString.append(_family+","+_genealogy+","+_generation+","+_birthday);
		return tempString.toString();
	}

	/* Lineage management ________________________________________________ */

	/**
	 * Called when creating an agent : update _generation and _genealogy field
	 */
	protected void recordGenealogy(Agent baby) {
		// Rob 18/1/11: Shuffled around slightly to include odd numbers
		baby._genealogy = _genealogy+ExtraMath.exp2(this._generation);

		this._generation++;
		baby._generation = this._generation;

		// Rob 25/1/11: we want to know if this happens
		if (baby._genealogy<0) {
			LogFile.writeLog("Warning: baby's genealogy has gone negative:");
			LogFile.writeLog("family "+baby._family+", genealogy "+baby._genealogy+", generation "+baby._generation);
		}

		// Rob 21/1/11: changed so that only the baby is given a new birthday
		// this._birthday = SimTimer.getCurrentTime();
		// baby._birthday = this._birthday;
		baby._birthday = SimTimer.getCurrentTime();
	}

	public String sendName(){
		return _family+"-"+_genealogy;
	}

	/**
	 * 
	 * 
	 */
	public void giveName() {
		_family = ++nextFamily;
	}

	public void setFamily(int family, int genealogy, int generation) {
		_family = family;
		_genealogy = genealogy;
		_generation = generation;
	}
}
