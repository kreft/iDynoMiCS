/**
 * Project iDynoMicS
 * ___________________________________________________________________________
 * Bacterium : 
 * ______________________________________________________
 * @since June 2006
 * @copyright -> see Idynomics.java
 * @version 0.1
 * @author Laurent Lardon
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;
import java.util.*;

import org.jdom.Element;

import idyno.SimTimer;
import simulator.agent.*;
import simulator.geometry.ContinuousVector;
import simulator.Simulator;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

public class MultiEpiBac extends BactEPS
{
	/* Parameters mutated from species parameters ___________________________ */
	/* Parameters specific to the agent _____________________________________ */

	public static StringBuffer plasmidTempString;

	// Plasmid hosted by this agent
	public Vector <MultiEpisome> plasmidHosted = new Vector<MultiEpisome>();

	//not being used in this version
	private Double _lastReception = 0.0;
	private Double _lastExchange = 0.0;
	private int _nCopy = 0;
	private int _status = 0;

	//sonia: 14-05-09
	//conjugation management
	public Vector <String>  partnerVector = new Vector <String>();
	public Vector <String>  plasmidVector = new Vector <String>();
	public static MultiEpiBac       _partner;

	protected LocatedAgent localAgent;

	//sonia: counting conjugation events
	//sonia: conjResult: only the donor cell will contain the information about the conjugation event: who was
	//its partner and where it is located
	public boolean conjResult = false;
	public static int conjugationEvents=0;

	//sonia 8-12-2010
	//distance based probability ordering management
	public Map<Double, LocatedAgent> test = new HashMap <Double, LocatedAgent> ();
	/* _________________________ CONSTRUCTOR _____________________________ */
	/**
	 * Empty constructor ; called to build a progenitor ; the speciesParameter
	 * have to be defined later
	 */
	public MultiEpiBac() {
		super();
		_speciesParam = new MultiEpiBacParam();
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		MultiEpiBac o = (MultiEpiBac) super.clone();
		o.plasmidHosted = new Vector <MultiEpisome>();
		MultiEpisome newEpisome;
		for (MultiEpisome anEpisome : plasmidHosted) {
			newEpisome = (MultiEpisome) anEpisome.clone();
			newEpisome.setHost(o);
			o.plasmidHosted.add(newEpisome);	
		}

		return o;
	}



	/**
	 * Called during species creation to build the progenitor
	 */
	@Override
	public void initFromProtocolFile(Simulator aSimulator, XMLParser aSpeciesRoot) {
		// Initialisation of the Located agent
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);

		// Create hosted plasmids
		for (String aSpeciesName : aSpeciesRoot.getChildrenNames("plasmid"))
			addPlasmid(aSpeciesName);

		// Genealogy and size management
		init();
	}

	//sonia 21/01/2011
	//this only works for the case where there's only one plasmid
	//however, could use the plasmidHosted vector to read in the names of the plasmids carried
	//everything else doesn't matter for my model (e.g. nCopy, lastReception or lastExchange)

	//public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
	// find the position to start at by using length and number of values read
	//int nValsRead = 3;
	//int iDataStart = singleAgentData.length - nValsRead;

	// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT

	// HGT parameters
	//_status        = Integer.parseInt(singleAgentData[iDataStart]);
	//_nCopy         = Integer.parseInt(singleAgentData[iDataStart+1]);
	//_lastReception = Double.parseDouble(singleAgentData[iDataStart+2]);
	//_lastExchange  = Double.parseDouble(singleAgentData[iDataStart+3]);

	// now go up the hierarchy with the rest of the data
	//String[] remainingSingleAgentData = new String[iDataStart];
	//for (int i=0; i<iDataStart; i++)
	//	remainingSingleAgentData[i] = singleAgentData[i];
	//super.initFromResultFile(aSim, remainingSingleAgentData);
	//}

	@Override
	public MultiEpiBac sendNewAgent() throws CloneNotSupportedException {
		MultiEpiBac baby = (MultiEpiBac) this.clone();
		baby.init();
		return baby;
	}


	@Override
	public void createNewAgent(ContinuousVector position) {
		try {
			// Clone the agent
			MultiEpiBac baby = sendNewAgent();
			baby.randomiseMass();
			baby.setLocation(position);
			baby.giveName();
			baby.updateSize();

			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			System.out.println("at createNewAgent in EpiBac error " + e);
		}
	}
	
	/* ______________________ CELL DIVISION ___________________ */

	@Override
	public void makeKid() throws CloneNotSupportedException 
	{
		/*
		 * Create the new instance.
		 */
		MultiEpiBac baby = sendNewAgent();
		/*
		 * Update the lineage.
		 */
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * Compute and apply movement to both cells.
		 */
		if ( ! Simulator.isChemostat )
		{
			setDivisionDirection(getInteractDistance(baby)/2);
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		/*
		 * For now the newborns' plasmids will be given the fixed copy number,
		 * that is they'll have the same copy number as their progenitor.
		 * Both daughters cells have cloned plasmids ; apply the segregation
		 */
		for ( int i = 0; i<plasmidHosted.size(); i++ )
			plasmidHosted.get(i).segregation(baby.plasmidHosted.get(i));
	}

	/* _____________________________ STEP __________________________________ */
	
	/**
	 * Method called by the STEP method (cf. the Agent class)
	 */
	@Override
	public void internalStep()
	{
		/*
		 * Check if some plasmid has a null copy number and remove it if
		 * necessary.
		 */
		checkMissingPlasmid();
		/*
		 * Everything else is the same as in Bacterium.
		 */
		super.internalStep();
	}


	public double fitnessCost (double initialCost, double rateOfDecay, double basalCost, double timeSpentInHost){

		double fCost = initialCost *(Math.exp(-(rateOfDecay*timeSpentInHost)))+ basalCost;

		return fCost;
	}
	
	/**
	 * TODO Is this any different from super?
	 */
	@Override
	public void grow()
	{
		Double deltaMass;

		int reacIndex;
		_netGrowthRate = 0.0;
		_netVolumeRate = 0.0;
		// Compute mass growth rate of each active reaction
		for (int iReac = 0; iReac<reactionActive.size(); iReac++)
		{
			// Compute the growth rate
			reacIndex = reactionActive.get(iReac);
			growthRate[reacIndex] = allReactions[reacIndex].computeMassGrowthRate(this);
			// Apply the growth rate on the particles
			for (int i = 0; i<particleYield[reacIndex].length; i++)
			{
				deltaMass = particleYield[reacIndex][i]*growthRate[reacIndex];
				_netGrowthRate += deltaMass;
				_netVolumeRate += deltaMass/getSpeciesParam().particleDensity[i];
				particleMass[i] += (deltaMass*SimTimer.getCurrentTimeStep());	
			}
		}
	}
	
	/**
	 * Remove a plasmid if its copy number reaches 0
	 */
	protected void checkMissingPlasmid()
	{
		Iterator<MultiEpisome> iter = plasmidHosted.iterator();
		while (iter.hasNext())
		{
			MultiEpisome anEpisome = iter.next();
			if ( anEpisome.getCopyNumber() <= 0 )
			{
				losePlasmid(anEpisome);
				anEpisome.die();
				iter.remove();
			}
		}
	}
	
	/**
	 * test if this cell can be a recipient for a given plasmid
	 * @param aPlasmid
	 * @return true if this cell is compatible
	 */
	public boolean isCompatible(MultiEpisome aPlasmid, MultiEpiBac partner) {

		boolean hMarker = false;
		boolean pMarker = false;

		//sonia:
		//if it is host specific, then we have to assess whether the plasmid to transfer (aPlasmid) can be
		// maintained in the recipient cell being queued (using the markers)	
		int sizeM = aPlasmid.getSpeciesParam().hostCompatibilityMarkers.size();		
		loopA:
			for (int i=0; i<sizeM; i++){
				if( partner.getName().equals(aPlasmid.getSpeciesParam().hostCompatibilityMarkers.get(i))){
					hMarker = true;
					break loopA;
				}
			}

		// sonia: 02.03.2010 test whether the partner contains any plasmids and if these are 
		// compatible (i.e., if they can replicate in the same cell) with the incoming plasmid

		int plListSize = partner.plasmidHosted.size();
		ArrayList<String> plHostedNames = new ArrayList<String>();

		if(plListSize>0){
			for (int i=0; i<plListSize; i++){
				plHostedNames.add(partner.plasmidHosted.get(i).getName());
			}


			int sizeP = aPlasmid.getSpeciesParam().plasmidCompatibilityMarkers.size();
			loopB:
				for (int i=0; i<sizeP; i++){
					for (int j=0; j< plHostedNames.size(); j++){
						if( plHostedNames.get(j).equals(aPlasmid.getSpeciesParam().plasmidCompatibilityMarkers.get(i))){
							pMarker = true;
							break loopB;
						}
					}
				}

		}else{
			pMarker = true;
		}

		if (hMarker & pMarker){
			return true;
		}else{
			return false;
		}


	}


	/**
	 * Search a recipient in your neighbourhood, and try to initiate a conjugation with it.
	 * 
	 */
	public void conjugate(double elapsedHGTtime) {

		// For each plasmid ready to conjugate search a number of potential recipients (partners) and conjugate

		//Randomise list of plasmids, specially useful in the incompatible plasmids scenario
		Collections.shuffle(plasmidHosted, ExtraMath.random);

		for (MultiEpisome aPlasmid : plasmidHosted) {

			//if(Simulator.isChemostat){

			/*		if(aPlasmid._newT == SimTimer.getCurrentIter()){
					//if this is a newly formed transconjugant, do nothing.
					// update its _newT field for the next round of conjugation that'll take place
					//in the next time step
					aPlasmid._newT = SimTimer.getCurrentIter()-1;

				}else{*/
			//searchConjugation(aPlasmid);


			//Biofilm
			//}else{
				
			
				//if (aPlasmid.isReadyToConjugate(elapsedHGTtime)){
					searchConjugation(aPlasmid);
			//}

			//}
		}

	}

	/* __________________ CONJUGATION ___________________________ */

	/**
	 * sonia: 01-05-09
	 * sonia: modified 8-12-2010
	 * 
	 * Run a dice to know if we initiate the conjugation (probability of transfer).
	 */

	public boolean acceptConjugation(MultiEpisome aPlasmid, MultiEpiBac partner, double distBasedProb) {

		double tP, rP;
		//sonia 8-12-2010
		boolean conjugate = true;

		/*		int plListSize = partner._plasmidHosted.size();
		ArrayList<String> plHostedNames = new ArrayList<String>();

		if(plListSize>0){
			for (int i=0; i<plListSize; i++){
			plHostedNames.add(partner._plasmidHosted.get(i).getName());
			}

		if (plListSize>0) {
			if(plHostedNames.contains(aPlasmid.getName())){
					conjugate = false;
				}else if (plHostedNames.contains("NHRa")){
					conjugate = false;
				}else if (plHostedNames.contains("NHRb")){
					conjugate = false;
				}else{
					conjugate = true;
				}
			}
		}else{
			conjugate = true;
		}*/

		tP = aPlasmid.getSpeciesParam().transferProb;

		//sonia: the recipient probability encompasses the retroransfer probability and enzyme restriction systems acting on the
		//recipient cell
		rP = partner.getSpeciesParam().recipientProbability;

		conjugate &= (ExtraMath.getUniRandDbl()<=tP*rP*distBasedProb); 

		return (conjugate);
	}



	public boolean testDonorTransfer(MultiEpisome aPlasmid){

		double tP;
		boolean conjugate = true;

		tP = aPlasmid.getSpeciesParam().transferProb;

		conjugate &= (ExtraMath.getUniRandDbl()<=tP); 

		return (conjugate);
	}



	/**
	 * \brief Initiate the search for a recipient cell in the neighbourhood.
	 * 
	 * @param aPlasmid
	 */
	public synchronized void searchConjugation(MultiEpisome aPlasmid)
	{
		if ( Simulator.isChemostat )
		{
			int i = ExtraMath.getUniRandInt(_agentGrid.agentList.size());
			SpecialisedAgent anAgent = _agentGrid.agentList.get(i);			
			if ( anAgent != this && anAgent instanceof MultiEpiBac)
			{
				_partner = (MultiEpiBac) anAgent;
				if ( isCompatible(aPlasmid, _partner) )
					acceptConjugation(aPlasmid, _partner, 1);
			}
		}
		else
		{
			/*
			 * Build your neighbourhood. If it's empty, nothing more to do.
			 */
			if ( aPlasmid.nbhList.isEmpty() ) 
				buildNbh(aPlasmid.getPilusLength(), aPlasmid);
			if ( aPlasmid.nbhList.isEmpty() )
				return;
			/*
			 * First test whether the plasmid will be transferred and then
			 * proceed with the distance-based probability of transfer to a
			 * nearest recipient.
			 */
			if ( testDonorTransfer(aPlasmid) )
			{
				Double cumProbSum = 0.0;
				for ( LocatedAgent agent : aPlasmid.nbhList )
					cumProbSum += agent._distCumProb;
				Double normRand = ExtraMath.getUniRandDbl()*cumProbSum;
				/*
				 * Find a neighbour to try conjugation with.
				 */
				LocatedAgent aLoc = null;
				for (int i = 0; i< aPlasmid.nbhList.size(); i++)
				{
					aLoc =	aPlasmid.nbhList.get(i);
					if( aLoc._distCumProb < normRand )
					{
						aLoc = aPlasmid.nbhList.remove(i);
						break;
					}
				}
				if ( aLoc != this && aLoc instanceof MultiEpiBac )
				{
					_partner = (MultiEpiBac) aLoc;
					if ( isCompatible(aPlasmid, _partner) )
						acceptConjugation(aPlasmid, _partner, 1);
				}
			}
		}
	}
	
	/**
	 * List all cells in a given neighbourhood : at the end of the method, the field
	 * listNbh contains all locatedAgents located in the neighbourhood
	 * @param nbhRadius
	 */
	public void buildNbh(Double nbhRadius, MultiEpisome aPlasmid)
	{
		/*
		 * Distance between cell surfaces.
		 */
		Double dist = 0.0;
		/*
		 * Radii of donor and recipient.
		 */
		Double donorRadius, recipRadius;
		/*
		 * Probability of conjugation success.
		 */
		Double distProb = 0.0;
		/*
		 * Search for potential recipients.
		 */
		Double radius = Math.ceil(nbhRadius/_agentGrid.getResolution());
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
		/*
		 * Now remove any agents that are too far (apply circular perimeter).
		 */
		for ( LocatedAgent aLocAgent : _myNeighbors )
		{
			if ( aLocAgent == this )
				continue;
			/*
			 * The distance between two cells is measured from their surface
			 * and not from the center of their mass.
			 */
			donorRadius = this.getRadius(false);
			recipRadius = aLocAgent.getRadius(false);
			dist = this.getDistance(aLocAgent) - donorRadius - recipRadius;
			if ( dist < nbhRadius )
			{
				distProb = ExtraMath.sq(donorRadius/(donorRadius+dist));
				aLocAgent._distProb = distProb;
				test.put(distProb, aLocAgent);	
			}	
		}
		_myNeighbors.clear();
		/*
		 * Order the recipients according to their distance to the donor cell.
		 */
		for ( Double reach : test.keySet() )
			if ( reach < nbhRadius )
				aPlasmid.nbhList.addLast(test.get(reach));
		/*
		 * Calculate and apply the cumulative probabilities.
		 */
		Double cumulative = 0.0;
		for ( LocatedAgent aLoc : aPlasmid.nbhList )
		{
			cumulative += aLoc._distProb;
			aLoc._distCumProb = cumulative;
		}
	}

	/* ______________________ HIGH LEVEL METHOD ___________________________ */

	/**
	 * Add a new plasmid to the list of hosted plasmids ; based on the
	 * species name of the plasmid.
	 */
	public void addPlasmid(String plasmidName)
	{
		try {
			MultiEpisome aPlasmid = (MultiEpisome)
					_species.getSpecies(plasmidName).sendNewAgent();
			plasmidHosted.add(aPlasmid);
			aPlasmid.setHost(this);
			/*
			 * When the cells carrying a plasmid are created we must set the
			 * field "timeSpentInHost" of the plasmid as being the time the
			 * cell was born (the current time).
			 */
			aPlasmid.timeSpentInHost = SimTimer.getCurrentTime();
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "EpiBac.addPlasmid("+plasmidName+")");
		}
	}
	
	/**
	 * 
	 * @param aPlasmid
	 */
	public void losePlasmid(MultiEpisome aPlasmid)
	{
		for (int aReaction : aPlasmid.reactionActive)
			removeReaction(allReactions[aReaction]);
	}
	
	@Override
	public MultiEpiBacParam getSpeciesParam()
	{
		return (MultiEpiBacParam) _speciesParam;
	}

	/**
	 * Used to write povray files
	 * 
	 * @return color of the host if empty, color of the first hosted plasmid
	 * else
	 */
	@Override
	public Color getColor()
	{
		MultiEpiBacParam param = getSpeciesParam();
		if ( plasmidHosted.isEmpty() )
			return param.rColor;
		for ( MultiEpisome plasmid : plasmidHosted )
			if ( plasmid.isTransConjugant() )
				return param.tColor;
		return param.dColor;
	}
	
	@Override
	public StringBuffer sendHeader()
	{
		return super.sendHeader().append(",plasmid,copyNumber,transconjugant");
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
		for ( MultiEpisome anEpi : plasmidHosted )
		{	
			tempString.append(",");
			tempString.append(anEpi.getSpeciesParam().plasmidName + ",");
			tempString.append(anEpi.getCopyNumber()+ ",");
			/*
			 * Count cells that carry a certain type of plasmid; the copy
			 * number of the plasmid is irrelevant for the time being -->
			 * that's done in the matlab scripts analysing the agent_sum xml
			 * files.
			 */
			tempString.append(anEpi.isTransConjugant() ? "1" : "0");
		}
		plasmidVector.clear();
		partnerVector.clear();
		return tempString;
	}

	@Override
	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		MultiEpiBacParam param = getSpeciesParam();

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