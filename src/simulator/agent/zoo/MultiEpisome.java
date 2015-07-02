/**
 * Project MicoDyna
 * Class describing the dynamic of one kind of plasmid in a cell
 * 
 * @since June 2006
 * @version 0.1
 * @author Laurent Lardon, DTU Denmark
 */

package simulator.agent.zoo;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.LinkedList;

import org.jdom.Element;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.agent.*;
import simulator.reaction.Reaction;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

public class MultiEpisome extends InfoAgent
{
	/**
	 * 
	 */
	protected MultiEpiBac              _host;

	/* Division management ___________________________________ */
	
	/**
	 * 
	 */
	private int _nCopy;

	/* Conjugation management __________________________________ */
	
	/**
	 * Already coded in EpisomeParam
	 * 
	 * Doesn't seem to be used here, but is in Episome!
	 */
	public Double lastExchange = 0.0;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double lastReception = 0.0;
	
	/**
	 * TODO Rename timeEnteredHost?
	 */
	public Double timeSpentInHost;
	
	/**
	 * This just gets the parameter from MultiEpisomeParam
	 * Already coded in EpisomeParam
	 */
	private Double _pilusLength;
	
	/**
	 * Doesn't seem to be used
	 */
	public Double _newT = 0.0;
	
	/* Reaction management ___________________________________ */
	
	/**
	 * Doesn't seem to be used
	 */
	private boolean             _conjugationIsOff = false;
	
	/**
	 * Doesn't seem to be used
	 */
	public boolean              isHot             = false;
	
	/**
	 * Doesn't seem to be used
	 */
	public Reaction[]            allReactions;
	
	/**
	 * Doesn't seem to be used
	 */
	protected ArrayList<Integer> reactionActive;
	
	/**
	 * Doesn't seem to be used
	 */
	protected ArrayList<Integer> reactionKnown;
	
	/**
	 * Array containing list of potential neighbours to be screened during the
	 * HGT time step.
	 * 
	 * TODO In Brian's code, this is done by EpiBac rather than Episome.
	 */
	protected LinkedList<LocatedAgent> nbhList = new LinkedList<LocatedAgent>();
	
	/* ____________________ CONSTRUCTOR _______________________________ */

	public MultiEpisome()
	{
		super();
		_speciesParam = new MultiEpisomeParam();
	}

	//sonia 12.10.09
	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException
	{
		MultiEpisome o = (MultiEpisome) super.clone();
		o._host = this._host;
		o._speciesParam = _speciesParam;
		o.allReactions = this.allReactions.clone();
		o.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
		o.reactionKnown = (ArrayList<Integer>) this.reactionKnown.clone();
		//sonia 11.10.2010
		o.nbhList = (LinkedList<LocatedAgent>) this.nbhList.clone();

		return o;
	}

	/**
     * Attributes the plasmid to an host (used for conjugation events and for cell division)
     * @param anHost
     */
	public void setHost(MultiEpiBac anHost)
	{
		lastReception = SimTimer.getCurrentTime();
		//lastExchange = -1;
		this._host = anHost;
		//sonia: have to decide on this... should we set the maximum copy number upon plasmid reception?
		//or should we model its replication/copy number control ?
		//setDefaultCopyNumber();
	}
	

	


	/* ______________________ CREATION _____________________________ */

	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp)
	{
		// Initialisation of the Located agent
		// super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		// init();
		
		_nCopy = getSpeciesParam().nCopy;
		_pilusLength = getSpeciesParam().pilusLength;
		int reacIndex;
		
		allReactions = aSim.reactionList;
		reactionKnown = new ArrayList<Integer>();
		reactionActive = new ArrayList<Integer>();
		
		for (Element aReactionMarkUp : xmlMarkUp.getChildrenElements("reaction"))
		{
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			// Add the reaction to the list of known (and active) reactions
			reactionKnown.add(reacIndex);
			if (aReactionMarkUp.getAttributeValue("status").equals("active"))
				reactionActive.add(reacIndex);
		}
		
	}
	
	/**
	 * Doesn't seem to be used
	 */
	public MultiEpisome reset(MultiEpisome aPlasmid)
	{
		aPlasmid._generation = 0;
		aPlasmid._genealogy = BigInteger.ZERO;
		aPlasmid.lastExchange = -1.0;
		aPlasmid.lastReception = -1.0;
		return aPlasmid;
	}
	
	/**
     * Used to initialise any new agent (progenitor or daughter cell)
     * 
     * @see sendNewAgent()
     */
	public void init()
	{
		// Lineage management : this is a new agent, he has no known parents
		_generation = 0;
		_genealogy = BigInteger.ZERO;
		//sonia 01.03.2010 changed from -1 to 0
		lastExchange = 0.0;
		lastReception = SimTimer.getCurrentTime();
	}

	/**
     * 
     */
	@Override
	public MultiEpisome sendNewAgent() throws CloneNotSupportedException
	{
		MultiEpisome baby = (MultiEpisome) this.clone();
		baby.init();
		return baby;
	}

	/**
     * 
     */
	@Override
	public void createNewAgent()
	{
		try
		{
			// Clone the plasmid
			MultiEpisome baby = this.sendNewAgent();
			// Register the plasmid (species population)
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "MultiEpisome.createNewAgent()");
		}
	}

	@Override
	public void registerBirth()
	{
		_species.notifyBirth();
	}

	/* __________________________ CELL DIVISION ____________________________ */
	
	/**
	 * Doesn't seem to be used
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		/*
		 * Clone the plasmid.
		 */
		MultiEpisome baby = this.sendNewAgent();
		/*
		 * Register the plasmid (species population).
		 */
		baby.registerBirth();
	}

	/* ____________________________________________________________________ */

	@Override
	public void internalStep()
	{
		/*
		 * For now we will only have a fixed number of plasmid copies, hence
		 * there will be no plasmid replication.
		 */
		//replication();
	}
	
	/**
	 * 
	 */
	public void die()
	{
		_species.notifyDeath();
	}
	
	/**
	 * 
	 * @return
	 */
	public Double getPilusLength()
	{
		return _pilusLength;
	}
	
	/**
	 * Called during the division of an EpiBac ; apply the segregation of
	 * plasmids, modify the number for the plasmid calling the method and
	 * sends the number for the other one.
	 * 
	 * @return
	 */
	public void segregation(MultiEpisome aPlasmid)
	{
		if ( ExtraMath.getUniRandDbl() > getSpeciesParam().lossProbability )
			_nCopy = 1;
		else
			_nCopy = 0;
		aPlasmid._nCopy = 1;
	}
	
	/**
	 * 
	 */
	@Override
	public MultiEpisomeParam getSpeciesParam()
	{
		return (MultiEpisomeParam) _speciesParam;
	}
	
	/**
	 * 
	 * @return
	 */
	public int getCopyNumber()
	{
		return _nCopy;
	}
	
	/**
	 * Now we can distinguish between vertical and horizontal transmission of
	 * the plasmid if the time of its reception is higher than the time at
	 * which its host was born, then the plasmid is harboured by a
	 * transconjugant cell.
	 * 
	 * @return
	 */
	public boolean isTransConjugant()
	{
		return this.lastReception > this._host.getBirthday();
	}
	
	// TODO
	@Override
	public StringBuffer sendHeader()
	{
		return new StringBuffer(
					"locationX,locationY,locationZ,mass,radius,growthRate");
	}

	// TODO
	@Override
	public StringBuffer writeOutput()
	{
		//tempString.append("copyNumber,transconjugant;");
		return new StringBuffer(_species.speciesName);
	}	
}
