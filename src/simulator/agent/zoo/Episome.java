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

import java.math.BigInteger;
import java.util.ArrayList;

import org.jdom.Element;

import idyno.SimTimer;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;
import simulator.agent.*;
import simulator.Simulator;

public class Episome extends InfoAgent
{
	/**
	 * The cell currently hosting this Episome.
	 */
	protected EpiBac _host;
	
	/**
	 * The number of plasmid copies of this episome in the host.
	 */
	protected int _nCopy;

	/* Conjugation management __________________________________ */
	public Double lastExchange, lastReception;
	
	protected ArrayList<Integer> reactionActive;
	
	/* ____________________ CONSTRUCTOR _______________________________ */

	public Episome() {
		super();
		_speciesParam = new EpisomeParam();
	}

	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException
	{
		Episome o = (Episome) super.clone();
		o._host = this._host;
		o._speciesParam = _speciesParam;

		o.reactionActive = (ArrayList<Integer>) this.reactionActive.clone();
		o.lastExchange = this.lastExchange;
		o.lastReception = this.lastReception;

		return o;
	}

	/**
	 * Attributes the plasmid to an host (used for conjugation events)
	 * 
	 * @param anHost
	 */
	public void setHost(EpiBac anHost)
	{
		lastReception = SimTimer.getCurrentTime();
		lastExchange = -1.0;
		_host = anHost;
		setDefaultCopyNumber();
	}

	public EpiBac getHost()
	{
		return _host;
	}

	/* ______________________ CREATION _____________________________ */

	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp)
	{
		// Initilaisation of the Located agent
		// super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		// init();
		_nCopy = getSpeciesParam().nCopy;
		int reacIndex;
		
		reactionActive = new ArrayList<Integer>();

		for (Element aReactionMarkUp : xmlMarkUp.getChildrenElements("reaction"))
		{
			reacIndex = aSim.getReactionIndex(
								aReactionMarkUp.getAttributeValue("name"));
			/*
			 * Add the reaction to the list of active reactions.
			 */
			if (aReactionMarkUp.getAttributeValue("status").equals("active"))
				reactionActive.add(reacIndex);
		}
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
		_genealogy = BigInteger.ZERO;
		lastExchange = -1.0;
		lastReception = -1.0;
	}

	/**
	 * 
	 */
	@Override
	public Episome sendNewAgent() throws CloneNotSupportedException
	{
		Episome baby = (Episome) this.clone();
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
			Episome baby = this.sendNewAgent();
			// Register the plasmid (species population)
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "Episome.createNewAgent()");
		}
	}
	
	/**
	 * TODO Is this any different from the super class? 
	 */
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
		Episome baby = this.sendNewAgent();
		/*
		 * Register the plasmid (species population).
		 */
		baby.registerBirth();
	}
	
	/* _______________________________________________________________________ */

	/**
	 * Test if the episome can be transfered.
	 */
	public Boolean isReadyToConjugate()
	{
		/*
		 * First check if the plasmid has sufficient copies.
		 */
		if ( _nCopy < 1 )
			return false;
		/*
		 * Now check timings: cannot conjugate if given/received a plasmid too
		 * recently.
		 */
		EpisomeParam param = getSpeciesParam();
		Double triggerTime = Math.max(param.exchangeLag + lastExchange,
										param.receptionLag + lastReception);
		return triggerTime >= SimTimer.getCurrentTime();
	}
	
	
	public void die()
	{
		_species.notifyDeath();
	}
	
	/**
	 * Returns the length (in um) of this plasmid's pilus.
	 */
	public Double getPilusRange()
	{
		return getSpeciesParam().pilusLength;
	}

	/**
	 * You are doing a conjugation! Update your lag variables.
	 */
	public void updateConjugationTime(Episome baby)
	{
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
	
	@Override
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
	}

	public Boolean isCompatible(Episome aPlasmid)
	{
		return aPlasmid.getSpeciesParam().compatibilityMarker !=
								this.getSpeciesParam().compatibilityMarker;
	}
}
