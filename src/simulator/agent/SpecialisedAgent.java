/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 */

/**
 * 
 * @since April 2007
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * 
 */

package simulator.agent;

import java.io.FileWriter;
import java.io.IOException;
import simulator.AgentContainer;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.geometry.ContinuousVector;
import utils.LogFile;
import utils.XMLParser;

public abstract class SpecialisedAgent extends Agent implements HasSpecies, Cloneable {

	// Parameters common to all agents of this class

	// Parameters common (strict equality) to all agents of the same Species
	protected Species        _species;
	public int               speciesIndex;
	public boolean           isDead = false;
	protected SpeciesParam   _speciesParam;
	protected AgentContainer _agentGrid;

	//sonia 26.04.2010
	//reason for agent's death
	public String death;

	// Parameters mutated from species parameters

	// Parameters specific to the agent

	/* ________________ CONSTRUCTOR ___________________________ */
	public SpecialisedAgent() {
		// Call constructor of parent class
		super();
		_speciesParam = new SpeciesParam();
	}

	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) {
		try {
			super.initFromProtocolFile(aSim, aSpeciesRoot);
			_agentGrid = aSim.agentGrid;
		} catch (Exception e) {
			LogFile.writeLog("Creating "+this.getSpecies().speciesName);
			System.exit(-1);
		}
	}

	public Object clone() throws CloneNotSupportedException {
		SpecialisedAgent o = (SpecialisedAgent) super.clone();

		// Copy the references (superficial copy)
		o._species = this._species;
		o._speciesParam = this._speciesParam;
		return (Object) o;
	}

	/**
	 * Create a new agent with mutated parameters based on species default
	 * values
	 */
	public abstract void createNewAgent();

	public abstract SpecialisedAgent sendNewAgent() throws CloneNotSupportedException;

	public void mutateAgent() {
		// Mutate parameters inherited
		super.mutateAgent();
		// Now mutate your parameters
	}

	public void mutatePop() {
		// Mutate parameters inherited
		// Now mutate your parameters
	}

	public void registerBirth() {
		_agentGrid = _species.currentSimulator.agentGrid;
		_agentGrid.registerBirth(this);
		_species.notifyBirth();
	}

	public void die(boolean isStarving) {
		// If you are too small, you must die !
		// Decrease the population of your species
		_species.notifyDeath();
		isDead = true;
		_agentGrid.registerDeath(this);
	}

	public double move() {
		return 0;
	}

	/* ______________ ACCESSORS & MUTATORS _______________________________ */

	public void setSpeciesParam(SpeciesParam aSpeciesParam) {
		_speciesParam = aSpeciesParam;
	}

	public SpeciesParam getSpeciesParam() {
		return _speciesParam;
	}

	public Species getSpecies() {
		return _species;
	}

	public void setSpecies(Species aSpecies) {
		_species = aSpecies;
		speciesIndex = aSpecies.speciesIndex;

	}

	public String getName() {
		return _species.speciesName;
	}

	public boolean willDetach() {
		return false;
	}

	public int getGridIndex() {
		return 0;
	}

	public double interact(boolean MUTUAL, boolean pull, boolean seq, double gain) {
		return 0;
	}

	public ContinuousVector followPressure(SoluteGrid pressure) {
		return new ContinuousVector(0, 0, 0);
	}

	public boolean isMoving() {
		return false;
	}

	/**
	 * this writes a color definition to the passed-in file; meant for later use in macros
	 * (Note that this routine is put here and not in Species to allow derived agents 
	 * to use different colors for different states; EpiBac is one example, with different colors
	 * for donor, recipient, and transconjugant states.)
	 * 
	 * @param theFile
	 */
	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		fr.write("#declare "+_species.speciesName+" = color rgb < ");
		fr.write(((float) _species.color.getRed()) / 255.0 + " , ");
		fr.write(((float) _species.color.getGreen()) / 255.0 + " , ");
		fr.write(((float) _species.color.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}
}
