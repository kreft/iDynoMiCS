package simulator.agent.zoo;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.LocatedAgent;
import simulator.agent.SpecialisedAgent;
import simulator.agent.Species;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Bacterium class that can host a number of Plasmids of differing
 * species. 
 * 
 * <p>Since this extends BactEPS, instances of this class may also produce and
 * excrete EPS.</p>
 * 
 * <p>This class is an amalgamation of the <b>EpiBac</b> class, written by
 * Brian Merkey, and the <b>MultiEpiBac</b> class, written by Sonia 
 * Martins.</p>
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidBac extends BactEPS
{
	/**
	 * Plasmids hosted by this bacterium.
	 */
	private LinkedList<Plasmid> _plasmidHosted = new LinkedList<Plasmid>();
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidBac()
	{
		super();
		_speciesParam = new PlasmidBacParam();
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		PlasmidBac out = (PlasmidBac) super.clone();
		out._plasmidHosted = new LinkedList<Plasmid>();
		Plasmid newPlasmid;
		for (Plasmid aPlasmid : _plasmidHosted )
		{
			newPlasmid = (Plasmid) aPlasmid.clone();
			out._plasmidHosted.add(newPlasmid);
		}
		return out;
	}
	
	/**
	 * Called during species creation to build the progenitor.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSimulator,
													XMLParser aSpeciesRoot)
	{
		/*
		 * Initialisation of the BactEPS, and its superclasses.
		 */
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		/*
		 * Create hosted plasmids.
		 */
		for ( String aSpeciesName : aSpeciesRoot.getChildrenNames("plasmid") )
			this.initPlasmid(aSpeciesName);
		/*
		 * Genealogy and size management.
		 */
		init();
		/*
		 * Finally, grab all the plasmid species names for reporting.
		 */
		this.collectPlasmidSpeciesNames(aSimulator);
	}


	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
	{
		/*
		 * First, grab all the plasmid species names for reporting.
		 */
		this.collectPlasmidSpeciesNames(aSim);
		/*
		 * Find the position to start at by using length and number of values
		 * read.
		 */
		int nValsRead = 3*getSpeciesParam().potentialPlasmids.size();
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * Read in info from the result file IN THE SAME ORDER AS IT WAS
		 * OUTPUT. HGT parameters:
		 */
		double r, d;
		int nCopy, spCounter = 0;
		for ( String plasmidName : getPotentialPlasmidNames() )
		{
			nCopy = Integer.parseInt(singleAgentData[iDataStart+3*spCounter]);
			if ( nCopy <= 0 )
				continue;
			r = Integer.parseInt(singleAgentData[iDataStart+3*spCounter + 1]);
			d = Integer.parseInt(singleAgentData[iDataStart+3*spCounter + 2]);
			Plasmid aPlasmid = this.initPlasmid(plasmidName);
			aPlasmid.setDetails(nCopy, r, d);
		}
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for ( int i = 0; i < iDataStart; i++ )
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}
	
	@Override
	public PlasmidBac sendNewAgent() throws CloneNotSupportedException
	{
		PlasmidBac baby = (PlasmidBac) this.clone();
		baby.init();
		return baby;
	}
	
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		/*
		 * Create the new instance and update the lineage.
		 */
		PlasmidBac baby = sendNewAgent();
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * Compute movement to apply to both cells.
		 */
		setDivisionDirection(getInteractDistance(baby) / 2);
		baby._movement.subtract(_divisionDirection);
		_movement.add(_divisionDirection);
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		/*
		 * Both daughters cells have cloned plasmids; apply the segregation.
		 * 
		 * TODO [Rob 31July2015]: I've followed Brian's logic here, but it
		 * seems strange that only the mother's plasmids are at risk of loss.
		 */
		for ( Plasmid aPlasmid : this._plasmidHosted )
			aPlasmid.applySegregation();
	}
	
	/*************************************************************************
	 * BASIC GETTERS & SETTERS
	 ************************************************************************/
	
	@Override
	public PlasmidBacParam getSpeciesParam()
	{
		return (PlasmidBacParam) _speciesParam;
	}
	
	/**
	 * \brief Provides a list of all Plasmids hosted by this PlasmidBac host.
	 * 
	 * @return LinkedList of Plasmid objects hosted by this PlasmidBac.
	 */
	public LinkedList<Plasmid> getPlasmidsHosted()
	{
		return _plasmidHosted;
	}
	
	/*************************************************************************
	 * 
	 ************************************************************************/
	
	@Override
	public void internalStep()
	{
		/*
		 * Check if any plasmids have illegal copy numbers.
		 */
		checkMissingPlasmid();
		/*
		 * BactEPS internalStep methods.
		 */
		grow();
		updateSize();
		manageEPS();
		if ( willDivide() )
			divide();
		if ( willDie() )
			die(true);
		/*
		 * Now try conjugating.
		 */
		conjugate();
	}
	
	/**
	 * \brief Remove any plasmids whose presence is no longer legal.
	 */
	protected void checkMissingPlasmid()
	{
		int nIter = this._plasmidHosted.size();
		Plasmid aPlasmid;
		for(int i = 0; i < nIter; i++)
		{
			aPlasmid = this._plasmidHosted.removeFirst();
			if ( aPlasmid.getCopyNumber() <= 0 )
				this.killPlasmid(aPlasmid);
			else
				this._plasmidHosted.addLast(aPlasmid);				
		}
		/*
		 * If any plasmids have been lost, refresh the reactions encoded by
		 * remaining plasmids: we should only remove those reactions  uniquely
		 * provided by the plasmid that was lost.
		 */
		if ( this._plasmidHosted.size() < nIter )
			this.refreshPlasmidReactions();
	}
	
	/**
	 * \brief Initialise a Plasmid to be hosted by this PlasmidBac.
	 * 
	 * <p>Note that the new Plasmid will have default copy number for the
	 * species it belongs to.</p>
	 * 
	 * @param plasmidName Species name of the new Plasmid.
	 */
	private Plasmid initPlasmid(String plasmidName)
	{
		Plasmid aPlasmid = null;
		try
		{
			aPlasmid = (Plasmid) 
							_species.getSpecies(plasmidName).sendNewAgent();
			this.welcomePlasmid(aPlasmid);
			aPlasmid.registerBirth();
			
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "PlasmidBac.initPlasmid("+plasmidName+")");
		}
		return aPlasmid;
	}
	
	/**
	 * \brief Receive a Plasmid into this host.
	 * 
	 * <p><b>[Rob 31July2015]</b> Removed updates to conjugation time, etc:
	 * This is now handled by the donor plasmid.</p>
	 * 
	 * @param aPlasmid Plasmid to be hosted by this PlasmidBac.
	 */
	public void welcomePlasmid(Plasmid aPlasmid)
	{
		this._plasmidHosted.add(aPlasmid);
		this.addPlasmidReactions(aPlasmid);
	}
	
	/**
	 * \brief Tell a Plasmid that it has been lost, and so should die.
	 * 
	 * <p><b>[Rob 31July2015]</b> Changed the part about losing reactions so
	 * that we now refresh all plasmid-encoded reactions. This avoids the
	 * possibility of a reaction being lost when another hosted plasmid still
	 * encodes for it.</p>  
	 * 
	 * @param plasmid Plasmid that was hosted by this PlasmidBac, but has now
	 * been lost.
	 */
	private void killPlasmid(Plasmid aPlasmid)
	{
		aPlasmid.die();
		LogFile.writeLog("Plasmid "+aPlasmid.sendName()+
											" lost from "+this.sendName());
	}
	
	/**
	 * \brief Add all the reactions encoded by a given Plasmid.
	 * 
	 * @param aPlasmid Plasmid whose reactions should be conferred to this
	 * PlasmidBac host.
	 */
	private void addPlasmidReactions(Plasmid aPlasmid)
	{
		for ( int reacIndex : aPlasmid.getReactionsEncoded() )
			this.addActiveReaction(allReactions[reacIndex], true);
	}
	
	/**
	 * \brief Lose all the reactions encoded by a given Plasmid.
	 * 
	 * <p><b>[Rob 31July2015]</b> Beware of losing reactions that are also
	 * encoded by other plasmids: don't assume that each encoded reaction is
	 * unique to a particular plasmid species!</p>
	 * 
	 * @param aPlasmid Plasmid whose reactions were conferred to this
	 * PlasmidBac host.
	 */
	private void losePlasmidReactions(Plasmid aPlasmid)
	{
		for ( int reacIndex : aPlasmid.getReactionsEncoded() )
			this.removeReaction(allReactions[reacIndex]);
	}
	
	/**
	 * \brief Refresh all the reactions conferred to this PlasmidBac host by
	 * its hosted Plasmids. 
	 */
	private void refreshPlasmidReactions()
	{
		for ( Plasmid aPlasmid : this._plasmidHosted )
			this.losePlasmidReactions(aPlasmid);
		for ( Plasmid aPlasmid : this._plasmidHosted )
			this.addPlasmidReactions(aPlasmid);
	}
	
	/**
	 * \brief Ask all your Plasmids to conjugate, if they are ready.
	 */
	protected void conjugate()
	{
		for ( Plasmid aPlasmid : this._plasmidHosted )
			if ( aPlasmid.isReadyToConjugate() )
				this.searchConjugation(aPlasmid);
	}
	
	/**
	 * List all cells in a given neighbourhood: at the end of the method, the
	 * field _myNeighbors contains all Bacteria with a surface-surface
	 * distance from this PlasmidBac of less than nbhRadius.
	 * 
	 * <p>Parameter <b>nbhRadius</b> is typically the pilus length.</p>
	 * 
	 * <p><b>[Rob 31July2016]</b> No need to shuffle this list: 
	 * LocatedAgent.pickNeighbour() uses a randomly generated integer for
	 * the index to return.</p>
	 * 
	 * @param nbhRadius double length (in um) of the maximum cell 
	 * surface-surface distance for another Bacterium to be considered a
	 * neighbor.
	 */
	public void buildNbh(double nbhRadius)
	{
		/*
		 * nbhRadius gives the distance OUTSIDE the donor agent that touches a
		 * recipient agent, and so we need to subtract the radii from
		 * getDistance() or add the radii to nbhRadius.
		 * 
		 * getDistance(aLocAgent) gets distance between cell centres.
		 */
		double temp = nbhRadius + this.getRadius(false);
		/*
		 * Find all neighbours in the Manhattan perimeter.
		 */
		this.getPotentialShovers(temp);
		/*
		 * Now remove agents that are too far (apply Euclidean perimeter).
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
	}
	
	/**
	 * \brief Search for partners and try to send them a plasmid.
	 * 
	 * @param aPlasmid A Plasmid, hosted by this PlasmidBac, that should try
	 * to conjugate with neighboring bacteria.
	 */
	public void searchConjugation(Plasmid aPlasmid)
	{
		/*
		 * Build a neighbourhood including only Bacteria. If this is empty,
		 * there is nothing more to do.
		 * 
		 * TODO [Rob 31July2015] In Sonia's code we only search once in the
		 * chemostat!
		 */
		if ( Simulator.isChemostat )
		{
			this._myNeighbors.clear();
			for ( SpecialisedAgent aSA : _agentGrid.agentList )
				if ( aSA instanceof Bacterium )
					this._myNeighbors.add((Bacterium) aSA);
		}
		else
			this.buildNbh( aPlasmid.getPilusRange() );
		if ( this._myNeighbors.isEmpty() )
			return;
		/*
		 * First update the plasmid's scan rate from its host's growth tone.
		 * The plasmid will calculate the number of neighbours it can
		 * look at (there may be some overflow from the previous timestep). 
		 */
		aPlasmid.updateScanRate(this.getScaledTone());
		/*
		 * Find a recipient(s) and try to send them a plasmid.
		 * 
		 * Note that LocatedAgent.pickNeighbor() picks a random member, with
		 * replacement, of this._myNeighbors.
		 */
		while ( aPlasmid.canScan() )
			aPlasmid.tryToSendPlasmid( this.pickNeighbor() );
	}
	
	/**
	 * \brief Growth tone as a linear interpolation between the two cutoffs
	 * specified in the protocol file.
	 * 
	 * <p>See Merkey <i>et al</i> (2011) p.5 for more details.</p>
	 * 
	 * @return double value in the range [0.0, 1.0]
	 */
	public double getScaledTone()
	{
		double out = getSpeciesParam().lowTonusCutoff; 
		out = (growthTone() - out)/(getSpeciesParam().highTonusCutoff - out);
		return Math.max(0.0, Math.min(1.0, out));
	}
	
	/**
	 * \brief Net growth rate as a fraction of the maximum rate.
	 * 
	 * <p>The value is as large as 1 but may also be negative.
	 * See Merkey <i>et al</i> (2011) p.5 for more details.</p>
	 * 
	 * TODO Is always using zero for the final index safe/correct here?  
	 * 
	 * @return Net growth rate divided by maximum growth rate.
	 */
	public Double growthTone()
	{
		Double tonusMax = 0.0;
		updateGrowthRates();
		for ( int reacIndex : reactionActive )
		{
			tonusMax += reactionKinetic[reacIndex][0] * 
												particleYield[reacIndex][0];
		}
		return _netGrowthRate / (tonusMax * particleMass[0]);
	}
	
	/*************************************************************************
	 * REPORTING
	 ************************************************************************/
	
	
	
	/**
	 * \brief Using the Simulator species list, collect the names of all
	 * Plasmid species that could be hosted by this PlasmidBac species.
	 * 
	 * @param aSim The Simulator this is running in.
	 */
	private void collectPlasmidSpeciesNames(Simulator aSim)
	{
		for ( Species aSpecies : aSim.speciesList )
		{
			if ( ! ( aSpecies.getProgenitor() instanceof Plasmid ) )
				continue;
			if ( ! ((Plasmid) aSpecies.getProgenitor()).isCompatible(this) )
				continue;
			getSpeciesParam().addPotentialPlasmidName(aSpecies.speciesName);
		}
	}
	
	/**
	 * \brief Provides a list of all plasmid species names that this
	 * PlasmidBac species could host.
	 * 
	 * @return ArrayList<String> of plasmid species names.
	 */
	private ArrayList<String> getPotentialPlasmidNames()
	{
		return this.getSpeciesParam().potentialPlasmids;
	}
	
	/**
	 * \brief Update the header for report output.
	 * 
	 * <p>For every plasmid species that could be hosted by this PlasmidBac
	 * species, the copy number, last reception time, and last donation time
	 * will always be reported, even if this is zero for some cells in this
	 * species.</p>
	 * 
	 * @see simulator.agent.LocatedAgent#sendHeader()
	 */
	@Override
	public StringBuffer sendHeader()
	{
		StringBuffer header = super.sendHeader();
		for (String plasmidSpeciesName :  getPotentialPlasmidNames())
		{
			header.append(","+plasmidSpeciesName+"CopyNumber");
			header.append(","+plasmidSpeciesName+"LastReception");
			header.append(","+plasmidSpeciesName+"LastDonation");
		}
		return header;
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
		int nCopy;
		double r, d;
		for (String plasmidSpeciesName :  getPotentialPlasmidNames() )
		{
			nCopy = 0;
			r = -Double.MAX_VALUE;
			d = -Double.MAX_VALUE;
			h: for ( Plasmid aPlasmid : _plasmidHosted )
				if ( aPlasmid.sendName().equals(plasmidSpeciesName) )
				{
					nCopy = aPlasmid.getCopyNumber();
					r = aPlasmid.getTimeRecieved();
					d = aPlasmid.getTimeLastDonated();
					break h;
				}
			tempString.append(","+nCopy+","+r+","+d);
		}
		return tempString;
	}
	
	/*************************************************************************
	 * POV-RAY
	 ************************************************************************/
	
	@Override
	public String getName()
	{
		StringBuffer out = new StringBuffer( _species.speciesName );
		//TODO
		return out.toString();
	}
	
	@Override
	public Color getColor()
	{
		//TODO
		return null;
	}
	
	@Override
	public void writePOVColorDefinition(FileWriter fr) throws IOException 
	{
		//TODO
	}
}
