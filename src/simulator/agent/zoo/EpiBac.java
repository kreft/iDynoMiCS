/**
 * Project iDynoMicS
 * ___________________________________________________________________________
 * @since June 2006
 * @copyright -> see Idynomics.java
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr)
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
 * ___________________________________________________________________________
 */

package simulator.agent.zoo;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import idyno.SimTimer;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;
import simulator.agent.*;
import simulator.Simulator;

// bvm 30.1.2009
// made this derive from BactEPS, not Bacterium
public class EpiBac extends BactEPS
{
	/* Parameters mutated from species parameters ___________________________ */

	/* Parameters specific to the agent _____________________________________ */
	// Plasmid hosted by this agent
	private LinkedList<Episome> _plasmidHosted = new LinkedList<Episome>();
	private Double _lastReception = 0.0;
	private Double _lastExchange = 0.0;
	private int _nCopy = 0;
	
	/**
	 * _status == 0: recipient
	 * _status == 1: donor
	 * _status == 11: transconjugant able to transfer
	 * _status == 101: donor has transferred
	 * _status == 111: transconjugant has transferred
	 */
	private int _status = 0;

	/* _________________________ CONSTRUCTOR _____________________________ */
	/**
	 * Empty constructor ; called to build a progenitor ; the speciesParameter
	 * have to be defined later
	 */
	public EpiBac()
	{
		super();
		_speciesParam = new EpiBacParam();
	}

	@Override
	public Object clone() throws CloneNotSupportedException
	{
		EpiBac out = (EpiBac) super.clone();
		out._plasmidHosted = new LinkedList<Episome>();
		Episome newEpisome;
		for (Episome anEpisome : _plasmidHosted)
		{
			newEpisome = (Episome) anEpisome.clone();
			newEpisome.setHost(out);
			out._plasmidHosted.add(newEpisome);
		}
		return (Object) out;
	}

	/**
	 * Called during species creation to build the progenitor.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSimulator,
													XMLParser aSpeciesRoot)
	{
		// Initialisation of the Located agent
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		// Create hosted plasmids
		for (String aSpeciesName : aSpeciesRoot.getChildrenNames("plasmid"))
			addPlasmid(aSpeciesName);
		// Genealogy and size management
		init();
	}


	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
	{
		/*
		 * Find the position to start at by using length and number of values
		 * read.
		 */
		int nValsRead = 4;
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * Read in info from the result file IN THE SAME ORDER AS IT WAS
		 * OUTPUT. HGT parameters:
		 */
		_status        = Integer.parseInt(singleAgentData[iDataStart]);
		_nCopy         = Integer.parseInt(singleAgentData[iDataStart+1]);
		_lastReception = Double.parseDouble(singleAgentData[iDataStart+2]);
		_lastExchange  = Double.parseDouble(singleAgentData[iDataStart+3]);
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for ( int i = 0; i < iDataStart; i++ )
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	@Override
	public EpiBac sendNewAgent() throws CloneNotSupportedException
	{
		EpiBac baby = (EpiBac) this.clone();
		baby.init();
		return baby;
	}

	/* ______________________ CELL DIVISION ___________________ */
	
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		// Create the new instance
		EpiBac baby = sendNewAgent();
		// Update the lineage
		recordGenealogy(baby);
		// Share mass of all compounds between two daughter cells and compute
		// new size
		divideCompounds(baby, getBabyMassFrac());
		// Compute movement to apply to both cells
		setDivisionDirection(getInteractDistance(baby) / 2);
		// move both daughter cells
		baby._movement.subtract(_divisionDirection);
		_movement.add(_divisionDirection);
		// Now register the agent inside the guilds and the agent grid
		baby.registerBirth();
		// Both daughters cells have cloned plasmids ; apply the segregation
		for (int i = 0; i < _plasmidHosted.size(); i++)
			_plasmidHosted.get(i).segregation(baby._plasmidHosted.get(i));
	}

	/* _____________________________ STEP __________________________________ */
	/**
	 * Method called by the STEP method (cf. the Agent class)
	 */
	@Override
	public void internalStep()
	{
		// Check if some plasmid has a null copy number and remove it if
		// necessary
		checkMissingPlasmid();
		// Compute mass growth over all compartments
		grow();
		updateSize();
		// test if the EPS capsule has to be excreted
		manageEPS();
		// Test division and death treshold
		if (willDivide())
			divide();
		if (willDie())
			die(true);
		// Try to conjugate
		try
		{
			conjugate();
			updateStatus();
		}
		catch(Exception e)
		{
			LogFile.writeError(e, "EpiBac.internalStep() conjugation step");
		}
	}

	/**
	 * Remove a plasmid if its copy number reaches 0.
	 */
	protected void checkMissingPlasmid()
	{
		int nIter = _plasmidHosted.size();
		Episome anEpi;
		for(int i = 0; i < nIter; i++)
		{
			anEpi = _plasmidHosted.removeFirst();
			if ( anEpi.getCopyNumber() <= 0 )
				losePlasmid(anEpi);
			else
				_plasmidHosted.addLast(anEpi);				
		}
	}

	/**
	 * Test if this cell can be a recipient for a given plasmid by checking
	 * whether any of the plasmids already hosted is incompatible with it.
	 * 
	 * @param aPlasmid
	 * @return true if this cell is compatible
	 */
	public Boolean isCompatible(Episome aPlasmid)
	{
		for ( Episome hostPlasmid : _plasmidHosted )
			if ( ! hostPlasmid.isCompatible(aPlasmid) )
				return false;
		return true;
	}

	/**
	 * Search a recipient in your neighbourhood, and try to initiate a
	 * conjugation with him.
	 * TODO : cells are attracted
	 */
	public void conjugate()
	{
		/*
		 * For each plasmid ready to conjugate, search a partner and conjugate.
		 */
		for ( Episome aPlasmid : _plasmidHosted )
			if ( aPlasmid.isReadyToConjugate() )
				searchConjugation(aPlasmid);
	}

	/* __________________ CONJUGATION ___________________________ */

	/**
	 * Search a partner and try to send him a plasmid.
	 * 
	 * @param aPlasmid
	 * @return issue of the transfer
	 */
	public void searchConjugation(Episome aPlasmid)
	{
		/*
		 * Build a neighbourhood including only Bacteria. If this is empty,
		 * there is nothing more to do.
		 */
		buildNbh( aPlasmid.getPilusRange() );
		if ( _myNeighbors.isEmpty() )
			return;
		/*
		 * Find a partner(s) and try to send them a plasmid.
		 */
		Double maxTest = getMaxTest();
		LocatedAgent aLoc;
		for ( int i = 0; i < maxTest; i++ )
		{
			aLoc = pickNeighbor();
			if ( aLoc instanceof EpiBac )
				tryToSendPlasmid(aPlasmid, (EpiBac) aLoc);
		}
	}

	/**
	 * List all cells in a given neighbourhood: at the end of the method, the
	 * field _myNeighbors contains all Bacteria with a surface-surface
	 * distance from this EpiBac of less than nbhRadius.
	 * 
	 * Parameter nbhRadius is typically the pilus length.
	 * 
	 * @param nbhRadius
	 */
	public void buildNbh(Double nbhRadius)
	{
		/*
		 * Manhattan perimeter.
		 */
		getPotentialShovers(nbhRadius);
		/*
		 * nbhRadius gives the distance OUTSIDE the donor agent that touches a
		 * recipient agent, and so we need to subtract the radii from
		 * getDistance() or add the radii to nbhRadius.
		 * 
		 * getDistance(aLocAgent) gets distance between cell centres.
		 */
		Double temp = nbhRadius + this.getRadius(false);
		/*
		 * Now remove too far agents (apply circular perimeter).
		 */
		LocatedAgent aLocAgent;
		for (int iter = 0; iter < _myNeighbors.size(); iter++)
		{
			aLocAgent = _myNeighbors.removeFirst();
			if ( aLocAgent == this || ! (aLocAgent instanceof Bacterium) )
				continue;
			if ( getDistance(aLocAgent) < (temp + aLocAgent.getRadius(false)) )
				_myNeighbors.addLast(aLocAgent);
		}
		Collections.shuffle(_myNeighbors, ExtraMath.random);
	}

	/**
	 * Test transfer proficiency and if success, send the plasmid.
	 * 
	 * @param aPlasmid
	 * @param aTarget
	 * @return issue of the test
	 */
	public void tryToSendPlasmid(Episome aPlasmid, EpiBac aTarget)
	{
		if ( aTarget.isCompatible(aPlasmid) && aPlasmid.testProficiency() )
		{
			EpiBac aHost = aPlasmid.getHost();
			String transType;
			int hostStatus = aHost.getStatus();
			if ( hostStatus == 1 )
				transType = "Donor initial transfer ";
			else if ( hostStatus == 101 )
				transType = "Donor subsequent transfer ";
			else if ( hostStatus == 11 )
				transType = "Transconjugant initial transfer ";
			else if ( hostStatus == 111 )
				transType = "Transconjugant subsequent transfer ";
			else
				transType= "Unknown transfer ("+hostStatus+") ";
			LogFile.writeLog(transType+aHost.sendName()+" -> "+aTarget.sendName()+
					" at time "+SimTimer.getCurrentTime());
			aTarget.receivePlasmid(aPlasmid);
		}
	}

	/* _______________________ HIGH LEVEL METHOD ____________________________ */
	
	/**
	 * Add a new plasmid to the list of hosted plasmids; based on the
	 * speciesName of the plasmid.
	 */
	public void addPlasmid(String plasmidName)
	{
		try
		{
			Episome aPlasmid = (Episome)
							_species.getSpecies(plasmidName).sendNewAgent();
			welcomePlasmid(aPlasmid);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "EpiBac.addPlasmid()");
		}
	}
	
	/**
	 * 
	 * @param aPlasmid
	 */
	public void receivePlasmid(Episome aPlasmid)
	{
		try
		{
			/*
			 * Create a new instance and set new plasmid descriptors (copy
			 * number, host, timers). Register the plasmid to the host.
			 */
			Episome baby = aPlasmid.sendNewAgent();
			welcomePlasmid(baby);
		}
		catch (Exception e)
		{
			utils.LogFile.writeLog("Error met in EpiBac.receivePlasmid()");
		}
	}
	
	public void welcomePlasmid(Episome aPlasmid)
	{
		aPlasmid.setHost(this);
		_plasmidHosted.add(aPlasmid);
		addPlasmidReaction(aPlasmid);
		aPlasmid.updateConjugationTime(aPlasmid);
	}
	
	/**
	 * 
	 * @param aPlasmid
	 */
	public void losePlasmid(Episome aPlasmid)
	{		
		for ( int aReaction : aPlasmid.reactionActive )
			removeReaction(allReactions[aReaction]);
		aPlasmid.die();
		LogFile.writeLog("Plasmid lost "+sendName());
	}

	/**
	 * Add active reaction coded on the plasmid to active reaction of the host.
	 * 
	 * @param aPlasmid
	 */
	public void addPlasmidReaction(Episome aPlasmid)
	{
		for ( int aReaction : aPlasmid.reactionActive )
			addActiveReaction(allReactions[aReaction], true);
	}

	
	@Override
	public EpiBacParam getSpeciesParam()
	{
		return (EpiBacParam) _speciesParam;
	}
	
	/**
	 * 
	 * 
	 * @return
	 */
	public Double getMaxTest()
	{
		Double lowTonus = getSpeciesParam().lowTonusCutoff; 
		Double highTonus = getSpeciesParam().highTonusCutoff;
		Double theTonus = sendTonus();
		Double scanSpeed = getSpeciesParam().scanSpeed;
		/*
		 * Too low, so return zero.
		 */
		if ( theTonus < lowTonus )
			scanSpeed = 0.0;
		/*
		 * Middle case, so do linear interpolation.
		 */
		else if ( theTonus < highTonus )
			scanSpeed *= (theTonus-lowTonus) / (highTonus-lowTonus);
		/*
		 * If neither of these is called we have a high tonus,
		 * so just return maximum (same effect as no growth dependence).
		 */
		// TODO Jan's notes suggested we should be using agentTimeStep instead
		// of the global timestep. Check this!
		//return scanSpeed * SimTimer.getCurrentTimeStep();
		return scanSpeed * _agentGrid.AGENTTIMESTEP;
	}


	/**
	 * Used to write povray files (replaces the version in LocatedAgent)
	 * 
	 * If this is...
	 * - a recipient, append "_r" to the species name.
	 * - a transconjugant, append "_t" to the species name.
	 * - a donot, append "_d" to the species name.
	 */
	@Override
	public String getName()
	{
		StringBuffer out = new StringBuffer( _species.speciesName );
		if ( _plasmidHosted.isEmpty() )
			out.append( "_r" );
		else if ( _plasmidHosted.getFirst().lastReception > _birthday )
			out.append( "_t" );
		else
			out.append( "_d" );
		return out.toString();
	}

	/**
	 * Used to write povray files
	 */
	@Override
	public Color getColor()
	{
		EpiBacParam param = getSpeciesParam();
		/*
		 * Recipients have no plasmid.
		 * Transconjugant received the plasmid after birth.
		 * Donor received the plasmid before/at birth.
		 */
		if ( _plasmidHosted.isEmpty() )
			return param.rColor;
		else if ( _plasmidHosted.getFirst().lastReception > _birthday )
			return param.tColor;
		else
			return param.dColor;
	}

	/* _______________ FILE OUTPUT _____________________ */

	/**
	 * Return the header file for this agent's values after sending those
	 * for super.
	 */
	@Override
	public StringBuffer sendHeader()
	{
		StringBuffer tempString = super.sendHeader();
		tempString.append(",status,copyNumber,lastReception,lastExchange");
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
	@Override
	public StringBuffer writeOutput()
	{
		StringBuffer tempString = super.writeOutput();
		// what about _plasmidHosted?
		tempString.append(","+_status+","+_nCopy+","+_lastReception+","+_lastExchange);
		return tempString;
	}
	
	/**
	 * This returns the net growth rate as a fraction of the maximum rate, 
	 * and the value is as large as 1 but may also be negative.
	 * 
	 * @return
	 */
	public Double sendTonus()
	{
		Double tonusMax = 0.0;
		for ( int reacIndex : reactionActive )
			tonusMax += reactionKinetic[reacIndex][0] * particleYield[reacIndex][0];
		return _netGrowthRate / (tonusMax * particleMass[0]);
	}
	
	/**
	 * 
	 */
	public void updateStatus()
	{
		_nCopy = 0;
		_status = 0;
		/*
		 * Unit digit specifies number of plasmids hosted.
		 * 
		 * TODO Rob 16June2015: potential problem here!
		 * If more than one plasmid is hosted and nCopy > 0, status will have
		 * a unit digit that is not 0 or 1. This will not be recognised by 
		 * tryToSendPlasmid()
		 * Not the end of the world, since it will only interfere with
		 * reporting, but still worth reconsidering.
		 */
		for ( Episome anEpi : _plasmidHosted )
		{
			_lastReception = Math.max(_lastReception, anEpi.lastReception);
			_lastExchange = Math.max(_lastExchange, anEpi.lastExchange);
			_nCopy = Math.max(_nCopy, anEpi.getCopyNumber());
			_status++;
		}
		if ( _nCopy == 0 )
			_status = 0;
		else
		{
			/*
			 * Transconjugant received the plasmid after birth;
			 * Donor received the plasmid before/at birth.
			 * Tens digit specifies if it is a transconjugant.
			 */
			if( _lastReception > _birthday )
				_status += 10;
			/*
			 * Hundreds digit specifies whether it has transferred a plasmid.
			 */
			if( _lastExchange > _birthday )
				_status += 100;
		}
	}
	
	/**
	 * Return the current status for use elsewhere.
	 */
	public int getStatus()
	{
		return _status;
	}

	/**
	 * Writes a color definition to the passed-in file; meant for later use in
	 * macros. Overrules the version in SpecialisedAgent.
	 * 
	 * @param theFile
	 */
	@Override
	public void writePOVColorDefinition(FileWriter fr) throws IOException
	{
		EpiBacParam param = getSpeciesParam();

		fr.write("#declare "+_species.speciesName+"_d = color rgb < ");
		fr.write((param.dColor.getRed()) / 255.0 + " , ");
		fr.write((param.dColor.getGreen()) / 255.0 + " , ");
		fr.write((param.dColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");

		fr.write("#declare "+_species.speciesName+"_r = color rgb < ");
		fr.write((param.rColor.getRed()) / 255.0 + " , ");
		fr.write((param.rColor.getGreen()) / 255.0 + " , ");
		fr.write((param.rColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");

		fr.write("#declare "+_species.speciesName+"_t = color rgb < ");
		fr.write((param.tColor.getRed()) / 255.0 + " , ");
		fr.write((param.tColor.getGreen()) / 255.0 + " , ");
		fr.write((param.tColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}
}