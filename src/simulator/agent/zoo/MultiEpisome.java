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

import utils.XMLParser;
import utils.ExtraMath;
import idyno.SimTimer;
import simulator.agent.*;
import simulator.reaction.Reaction;
import simulator.Simulator;

public class MultiEpisome extends InfoAgent
{
	//sonia:I've changed it to public
	public MultiEpiBac              _host;

	/* Division management ___________________________________ */
	
	//sonia:I've changed it to public
	private int _nCopy;

	/* Conjugation management __________________________________ */
	public double               lastExchange=0, lastReception=0, timeSpentInHost;
	private double              _pilusLength;
	//public boolean				_newT = false;
   //sonia 01.03.2010
	public double _newT =0;
	/* Reaction management ___________________________________ */
	private boolean             _conjugationIsOff = false;
	public boolean              isHot             = false;

	public Reaction[]            allReactions;
	protected ArrayList<Integer> reactionActive;
	protected ArrayList<Integer> reactionKnown;
	
	//sonia 11.10.2010 array containing list of potential nbh to be screened during the hgt time step
	protected LinkedList<LocatedAgent> nbhList = new LinkedList<LocatedAgent>();


	
	
	/* ____________________ CONSTRUCTOR _______________________________ */

	public MultiEpisome() {
		super();
		_speciesParam = new MultiEpisomeParam();
	}

	//sonia 12.10.09
	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException {
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
	public void setHost(MultiEpiBac anHost) {
		
		//sonia 4.10.2010 
	
		lastReception = SimTimer.getCurrentTime();
		//lastExchange = -1;
		this._host = anHost;
		//sonia: have to decide on this... should we set the maximum copy number upon plasmid reception?
		//or should we model its replication/copy number control ?
		//setDefaultCopyNumber();
	}
	

	


	/* ______________________ CREATION _____________________________ */

	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) {
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

	public MultiEpisome reset(MultiEpisome aPlasmid){
		aPlasmid._generation = 0;
		aPlasmid._genealogy = BigInteger.ZERO;
		aPlasmid.lastExchange = -1;
		aPlasmid.lastReception = -1;
		
		return aPlasmid;
	}
	
	/**
     * Used to initialize any new agent (progenitor or daughter cell)
     * @see sendNewAgent()
     */
	public void init() {
		// Lineage management : this is a new agent, he has no known parents
		_generation = 0;
		_genealogy = BigInteger.ZERO;
		//sonia 01.03.2010 changed from -1 to 0
		lastExchange = 0;
		lastReception = SimTimer.getCurrentTime();
	}

	/**
     * 
     */
	@Override
	public MultiEpisome sendNewAgent() throws CloneNotSupportedException {
		MultiEpisome baby = (MultiEpisome) this.clone();
		baby.init();
		return baby;
	}

	/**
     * 
     */
	@Override
	public void createNewAgent() {
		try {
			// Clone the plasmid
			MultiEpisome baby = this.sendNewAgent();
			// Register the plasmid (species population)
			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			System.out.println("at Episome: createNewAgent error " + e);
		}
	}

	@Override
	public void registerBirth()
	{
		_species.notifyBirth();
	}

	/* __________________________ CELL DIVISION ____________________________ */

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
