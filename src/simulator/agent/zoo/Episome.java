/**
 * Project iDynoMicS
 * ______________________________________________________
 * @since June 2006
 * @copyright -> see Idynomics.java
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr)
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;

import java.util.ArrayList;

import org.jdom.Element;

import idyno.SimTimer;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;
import simulator.agent.*;
import simulator.reaction.Reaction;
import simulator.Simulator;

public class Episome extends InfoAgent
{
	protected EpiBac _host;
	
	/* Division management ___________________________________ */
	protected int _nCopy;

	/* Conjugation management __________________________________ */
	public Double lastExchange, lastReception;
	protected Double _pilusLength;

	/* Reaction management ___________________________________ */
	protected Boolean _isRepressed = false;
	public Boolean isHot = false;

	public Reaction[] allReactions;
	protected ArrayList<Integer> reactionActive;
	protected ArrayList<Integer> reactionKnown;

	/* ____________________ CONSTRUCTOR _______________________________ */

	public Episome() {
		super();
		_speciesParam = new EpisomeParam();
	}

	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException {
		Episome o = (Episome) super.clone();
		o._host = this._host;
		o._speciesParam = _speciesParam;

		o.allReactions = this.allReactions.clone();
		o.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
		o.reactionKnown = (ArrayList<Integer>) this.reactionKnown.clone();
		o.lastExchange = this.lastExchange;
		o.lastReception = this.lastReception;

		return (Object) o;
	}

	/**
	 * Attributes the plasmid to an host (used for conjugation events)
	 * 
	 * @param anHost
	 */
	public void setHost(EpiBac anHost) {
		lastReception = SimTimer.getCurrentTime();
		lastExchange = -1.0;
		_host = anHost;
		setDefaultCopyNumber();
	}

	public EpiBac getHost() {
		return _host;
	}

	/* ______________________ CREATION _____________________________ */

	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) {
		// Initilaisation of the Located agent
		// super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		// init();
		_nCopy = getSpeciesParam().nCopy;
		_pilusLength = getSpeciesParam().pilusLength;
		int reacIndex;

		allReactions = aSim.reactionList;
		reactionKnown = new ArrayList<Integer>();
		reactionActive = new ArrayList<Integer>();

		for (Element aReactionMarkUp : xmlMarkUp.buildSetMarkUp("reaction")) {
			reacIndex = aSim.getReactionIndex(aReactionMarkUp
					.getAttributeValue("name"));

			// Add the reaction to the list of known (and active) reactions
			reactionKnown.add(reacIndex);
			if (aReactionMarkUp.getAttributeValue("status").equals("active")) {
				reactionActive.add(reacIndex);
			}
		}
	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// this writes no unique values, so doesn't need unique reading-in
		// (for a template on how to read in data, look in LocatedAgent.java)
		super.initFromResultFile(aSim,singleAgentData);
	}

	/**
	 * Used to initialize any new agent (progenitor or daughter cell)
	 * 
	 * @see sendNewAgent()
	 */
	public void init()
	{
		// Lineage management : this is a new agent, he has no known parents
		_generation = 0;
		_genealogy = 0;
		lastExchange = -1.0;
		lastReception = -1.0;
	}

	/**
	 * 
	 */
	public Episome sendNewAgent() throws CloneNotSupportedException
	{
		Episome baby = (Episome) this.clone();
		baby.init();
		return baby;
	}

	/**
	 * 
	 */
	public void createNewAgent()
	{
		try
		{
			// Clone the plamid
			Episome baby = this.sendNewAgent();
			// Mutate parameters
			baby.mutatePop();
			// Register the plasmid (species population)
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "Episome.createNewAgent()");
		}
	}

	public void mutatePop()
	{
		// Mutate inherited parameters
		super.mutatePop();
		// Now mutate your parameters
	}

	public void registerBirth()
	{
		_species.notifyBirth();
	}

	/* __________________________ CELL DIVISION ____________________________ */

	public void makeKid() throws CloneNotSupportedException
	{
		// Clone the plamid
		Episome baby = this.sendNewAgent();
		
		// Mutate parameters
		baby.mutateAgent();
		
		// Register the plasmid (species population)
		baby.registerBirth();
	}

	public void mutateAgent()
	{
		// Mutate inherited parameters
		super.mutateAgent();
		// Now mutate your parameters
	}

	/* _______________________________________________________________________ */

	public void internalStep()
	{
		
	}

	/**
	 * Test if the episome can be transfered.
	 */
	public boolean isReadyToConjugate()
	{
		EpisomeParam param = getSpeciesParam();
		isHot = false;
		
		// not enough copys
		if (_isRepressed)
			return false;
		if (_nCopy < 1)
			return false;
		
		// You have given a plasmid a few minutes ago
		if (SimTimer.getCurrentTime() - lastExchange < param.exchangeLag)
			return false;
		
		// You have received this plasmid a few minutes ago
		if (SimTimer.getCurrentTime() - lastReception < param.receptionLag)
			return false;
		
		return true;
	}

	public void die()
	{
		_species.notifyDeath();
	}

	public double getPilusRange()
	{
		return _pilusLength;
	}

	/**
	 * You are doing a conjugation ! Update your parameters
	 */
	public void givePlasmid(Episome baby)
	{
		// Update your time counter for conjugation-able recovery
		lastExchange = SimTimer.getCurrentTime();
		baby.lastReception = SimTimer.getCurrentTime();		
	}

	/**
	 * Called during the division of an EpiBac ; apply the segregation of
	 * plasmids, modify the number for the plasmid calling the method and sends
	 * the number for the other one.
	 * 
	 * @param aPlasmid Episome 
	 */
	public void segregation(Episome aPlasmid)
	{
		if (ExtraMath.getUniRandDbl() > getSpeciesParam().lossProbability)
		{
			_nCopy = 1;
			aPlasmid._nCopy = 1;
		}
		else
		{
			_nCopy = 0;
			aPlasmid._nCopy = 1;
		}
	}

	public EpisomeParam getSpeciesParam()
	{
		return (EpisomeParam) _speciesParam;
	}

	public int getCopyNumber()
	{
		return _nCopy;
	}

	public void setDefaultCopyNumber()
	{
		_nCopy = getSpeciesParam().nCopy;
	}

	public Boolean testProficiency()
	{
		Double alea = ExtraMath.getUniRandDbl();
		return (alea <= getSpeciesParam().transferProficiency);
		
		// previous (LAL) growth-dependence mechanism
		//double prof = getSpeciesParam().transferProficiency;
		//if (_host.sendTonus() < 0.25)
		//	prof = prof / 5;
		//return (alea <= prof);
	}

	public Boolean isCompatible(Episome aPlasmid)
	{
		return aPlasmid.getSpeciesParam().compatibilityMarker != this
				.getSpeciesParam().compatibilityMarker;
	}
	
	public int giveStatus()
	{
		return 1;
	}
	
	/* _______________ FILE OUTPUT _____________________ */
	
	public StringBuffer sendHeader()
	{
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = super.sendHeader();
		//tempString.append(",");
		// _host, _isRepressed, _isHot, lastExchange, lastReception
		return tempString;
	}
	
	/**
	 * \brief Creates an output string of information generated on this
	 * particular agent.
	 * 
	 * Used in creation of results files.
	 * Writes the data matching the header file.
	 * 
	 * @return	String containing results associated with this agent.
	 */
	public StringBuffer writeOutput()
	{
		StringBuffer tempString = super.writeOutput();
		return tempString;
	}
}
