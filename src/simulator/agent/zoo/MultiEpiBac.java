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
import java.util.*;

import idyno.SimTimer;

import org.jdom.Element;

import utils.ExtraMath;
import utils.XMLParser;


import simulator.Simulator;
import simulator.agent.*;
import simulator.geometry.ContinuousVector;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;

public class MultiEpiBac extends BactEPS {

	/* Parameters mutated from species parameters ___________________________ */
	/* Parameters specific to the agent _____________________________________ */

	public static StringBuffer plasmidTempString;

	// Plasmid hosted by this agent
	public Vector <MultiEpisome> _plasmidHosted   = new Vector <MultiEpisome>();

	//not being used in this version
	private double _lastReception=0;
	private double _lastExchange=0;
	private int _nCopy=0;
	private int _status=0;

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
	public Map<Double, LocatedAgent> teste = new HashMap <Double, LocatedAgent> ();
	/* _________________________ CONSTRUCTOR _____________________________ */
	/**
	 * Empty constructor ; called to build a progenitor ; the speciesParameter
	 * have to be defined later
	 */
	public MultiEpiBac() {
		super();
		_speciesParam = new MultiEpiBacParam();
	}

	public Object clone() throws CloneNotSupportedException {
		MultiEpiBac o = (MultiEpiBac) super.clone();
		o._plasmidHosted = new Vector <MultiEpisome>();
		MultiEpisome newEpisome;
		for (MultiEpisome anEpisome : _plasmidHosted) {
			newEpisome = (MultiEpisome) anEpisome.clone();
			newEpisome.setHost(o);
			o._plasmidHosted.add(newEpisome);	
		}

		return (Object) o;
	}



	/**
	 * Called during species creation to build the progenitor
	 */
	public void initFromProtocolFile(Simulator aSimulator, XMLParser aSpeciesRoot) {
		// Initialisation of the Located agent
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);

		// Create hosted plasmids
		for (Element aSpeciesMarkUp : aSpeciesRoot.buildSetMarkUp("plasmid")) {
			addPlasmid(aSpeciesMarkUp.getAttributeValue("name"));		
		}

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

	public MultiEpiBac sendNewAgent() throws CloneNotSupportedException {
		MultiEpiBac baby = (MultiEpiBac) this.clone();
		baby.init();
		return baby;
	}


	public void createNewAgent(ContinuousVector position) {
		try {
			// Clone the agent
			MultiEpiBac baby = sendNewAgent();
			baby.mutatePop();
			baby.setLocation(position);
			baby.giveName();
			baby.updateSize();

			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			System.out.println("at createNewAgent in EpiBac error " + e);
		}
	}

	public void mutatePop() {
		// Mutate inherited parameters
		super.mutatePop();
		// Now mutate your parameters
	}

	/* ______________________ CELL DIVISION ___________________ */

	public void mutateAgent() {
		// Mutate inherited parameters
		super.mutateAgent();


		// Now mutate your parameters
	}

	public void makeKid() throws CloneNotSupportedException {
		// Create the new instance
		MultiEpiBac baby = (MultiEpiBac) sendNewAgent();
		baby.mutateAgent();

		// Update the lineage
		recordGenealogy(baby);

		// Share mass of all compounds between two daughter cells and compute
		// new size
		divideCompounds(baby, getBabyMassFrac());

		if(Simulator.isChemostat){

		}else{
			// Compute movement to apply to both cells
			setDivisionDirection(getInteractDistance(baby)/2);

			// move both daughter cells
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}

		// Now register the agent inside the guilds and the agent grid
		baby.registerBirth();

		//sonia 6/10/09
		//for now the newborns' plasmids will be given the fixed copy number,
		// that is they'll have the same copy number as their progenitors

		// Both daughters cells have cloned plasmids ; apply the segregation
		for (int i = 0; i<_plasmidHosted.size(); i++) {
			_plasmidHosted.get(i).segregation(baby._plasmidHosted.get(i));

		}

	}

	/* _____________________________ STEP __________________________________ */
	/**
	 * Method called by the STEP method (cf. the Agent class)
	 */
	public void internalStep() {

		// Check if some plasmid has a null copy number and remove it if
		// necessary
		checkMissingPlasmid();

		// Compute mass growth over all compartments
		grow();


		//sonia 11.10.2010 the hgt will be carried out in a separate function

		/*conjugate();

		conjResult = false;
		plasmidVector.clear();
		partnerVector.clear();*/

		// test if the EPS capsule has to be excreted
		updateSize();	

		manageEPS();

		// Apply this mass growth of all compounds on global radius and mass

		// Divide if you have to
		if (willDivide()) divide();

		// Die if you have to

		if (willDie()) die(true);


	}


	public double fitnessCost (double initialCost, double rateOfDecay, double basalCost, double timeSpentInHost){

		double fCost = initialCost *(Math.exp(-(rateOfDecay*timeSpentInHost)))+ basalCost;

		return fCost;
	}




	public void grow() {
		double deltaMass;
		int reacIndex;
		_netGrowthRate = 0;
		_netVolumeRate = 0;
		// Compute mass growth rate of each active reaction
		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			// Compute the growth rate
			reacIndex = reactionActive.get(iReac);
			growthRate[reacIndex] = allReactions[reacIndex].computeMassGrowthRate(this);

			// Apply the growth rate on the particles
			for (int i = 0; i<particleYield[reacIndex].length; i++) {
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
	protected void checkMissingPlasmid() {

		Iterator<MultiEpisome> iter = _plasmidHosted.iterator();

		while (iter.hasNext()) {
			MultiEpisome anEpisome = iter.next();
			if (anEpisome.getCopyNumber()<=0) {
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

		int plListSize = partner._plasmidHosted.size();
		ArrayList<String> plHostedNames = new ArrayList<String>();

		if(plListSize>0){
			for (int i=0; i<plListSize; i++){
				plHostedNames.add(partner._plasmidHosted.get(i).getName());
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
		Collections.shuffle(_plasmidHosted, ExtraMath.random);

		for (MultiEpisome aPlasmid : _plasmidHosted) {

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

	/**
	 * Remove agent and all references from the system
	 */
	public void die(boolean isStarving) {
		super.die(isStarving);
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
	 * Initiate the search for a recipient cell in the neighbourhood.
	 * @param aPlasmid
	 */
	public synchronized boolean searchConjugation(MultiEpisome aPlasmid) {

		//	int cellExamined = 0;
		LocatedAgent aLoc = null;
		// Search a compatible recipient in your neighbourhood
		Boolean test = false;


		//sonia 21.10.09
		// the scan speed should be multiplied by the agent time step, otherwise the cell is
		// screening more agents than  it should per unit of time; if the time step is higher than 1h,
		// then we also have to account for that by adjusting the number of screened individuals
		// per agentTimestep per hour of the global time step

		/*		double maxTest = 0;

		if(SimTimer.getCurrentTimeStep()>1){
			maxTest = (this.getSpeciesParam().scanSpeed * _agentGrid.AGENTTIMESTEP)/SimTimer.getCurrentTimeStep();
		}else{
			maxTest = this.getSpeciesParam().scanSpeed * _agentGrid.AGENTTIMESTEP;
		}*/

		if(Simulator.isChemostat){

			Collections.shuffle(_agentGrid.agentList, ExtraMath.random);
				Agent anAgent = _agentGrid.agentList.getFirst();
			
				//for(Agent anAgent: _agentGrid.agentList){
				//if (aPlasmid.isReadyToConjugate()){ //sonia 11.10.2010 no need to check this because
													  // only one agent can be infected during one HGT step
							
					//if(anAgent instanceof Bacterium){ 
					//if (cellExamined>maxTest) break;
				
						if (anAgent instanceof MultiEpiBac) {
							if (anAgent != this){
								_partner = (MultiEpiBac) anAgent;
								test = isCompatible(aPlasmid, _partner);
								test &= acceptConjugation(aPlasmid, _partner,1);
/*								if (test){
								sendPlasmid(aPlasmid, _partner, elapsedHGTtime);	
								}*/
							}
						}
					
			return false;

		}else{ //Biofilm

/*			if(elapsedHGTtime == SimTimer.getCurrentTime()){
				aPlasmid.nbhList.clear();
			}*/


			if (aPlasmid.nbhList.isEmpty()){
				// Build your neighbourhood
				buildNbh(aPlasmid.getPilusLength(), aPlasmid);
			}


			//if(_myNeighbors.isEmpty()) return false;
			if(aPlasmid.nbhList.isEmpty()) return false;


			//sonia 8-12-2010
			// First test whether the plasmid will be transferred and then proceed with the 
			// distance-based probability of transfer to a nearest recipient

			if(testDonorTransfer(aPlasmid)){


				double cumProbSum=0;

				for (int i=0; i< aPlasmid.nbhList.size(); i++){

					aLoc =	aPlasmid.nbhList.get(i);
					cumProbSum += aLoc._distCumProb;
				}

				double random = 0;
				double normRand =0;
				random = ExtraMath.getUniRandDbl();
				normRand = random*cumProbSum;

				int pos=0;
				for (int i=0; i< aPlasmid.nbhList.size(); i++){

					aLoc =	aPlasmid.nbhList.get(i);

					if(aLoc._distCumProb<normRand){
						aLoc = aPlasmid.nbhList.remove(i);
						pos=i;
						break;
					}
				}


				//if (aPlasmid.isReadyToConjugate()){				
				//	if (aLoc instanceof Bacterium) cellExamined++;
				//	if (cellExamined>maxTest) break;

				// aLoc = aPlasmid.nbhList.remove(pos);
				if (aLoc instanceof MultiEpiBac) {
					if(aLoc != this){
						_partner = (MultiEpiBac) aLoc;
						test = isCompatible(aPlasmid, _partner);
						test &= acceptConjugation(aPlasmid, _partner, normRand);

//						if (test) sendPlasmid(aPlasmid, _partner, elapsedHGTtime);
					}

				}
			}

			return false;
		}



	}

	/**
	 * List all cells in a given neighbourhood : at the end of the method, the field
	 * listNbh contains all locatedAgents located in the neighbourhood
	 * @param nbhRadius
	 */
	public void buildNbh(double nbhRadius, MultiEpisome aPlasmid) {

		double dist=0;
		double dRadius=0;
		double rRadius=0;		
		double distProb=0;

		int radius = (int) Math.ceil(nbhRadius/_agentGrid.getResolution());
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);


		// Now remove too far agents (apply circular perimeter)
		for (int iter = 0; iter<_myNeighbors.size(); iter++) {

			LocatedAgent aLocAgent = _myNeighbors.removeFirst();
			dist = this.getDistance(aLocAgent);

			//sonia 4/10/2010
			// the distance between two cells is measured from their surface and not from the centre of their mass
			dRadius=this.getRadius(false);
			rRadius = aLocAgent.getRadius(false);
			dist = dist - dRadius - rRadius;

			if(dist<nbhRadius){

				aLocAgent._distProb = (dRadius*dRadius)/ ((dRadius + dist) * (dRadius + dist));
				distProb = aLocAgent._distProb ;

				//aPlasmid.nbhList.addLast(aLocAgent);
				//_myNeighbors.addLast(aLocAgent);
				teste.put(distProb, aLocAgent);	
			}	
		}

		//sonia 10.2010 code for ordering the recipients according to their distance to the donor cell

		for (Iterator<Double> iter1 = teste.keySet().iterator(); iter1.hasNext();) {

			double reach = iter1.next();
			//System.out.println("porbability based distance is " + reach);
			if (reach<nbhRadius){
				aPlasmid.nbhList.addLast(teste.get(reach));
			}
		}
		//System.out.println("nbhList with close enough recipients size is " + aPlasmid.nbhList.size());

		double previousVal=0;
		double newVal=0;
		int counter =1;

		for (int i = 0; i< aPlasmid.nbhList.size(); i++){

			if (counter ==1){
				LocatedAgent aLocAgentNew = aPlasmid.nbhList.get(i);

				newVal = aLocAgentNew._distProb;

				aLocAgentNew._distCumProb = newVal;

			}else{

				LocatedAgent aLocAgentPrev = aPlasmid.nbhList.get(i-1);
				LocatedAgent aLocAgentNew = aPlasmid.nbhList.get(i);

				previousVal = aLocAgentPrev._distProb;
				newVal = aLocAgentNew._distProb;

				aLocAgentNew._distCumProb = previousVal + newVal;
			}
			counter++;

		}

		//Collections.shuffle(_myNeighbors);
		//Collections.shuffle(aPlasmid.nbhList);

	}

	/* _______________________ HIGH LEVEL METHOD ____________________________ */


	/**
	 * Add a new plasmid to the list of hosted plasmids ; based on the
	 * species name of the plasmid.
	 */
	public void addPlasmid(String plasmidName) {
		try {
			MultiEpisome aPlasmid = (MultiEpisome) _species.getSpecies(plasmidName).sendNewAgent();
			_plasmidHosted.add(aPlasmid);
			aPlasmid.setHost(this);
			//sonia: 08-06-09
			//when the cells carrying a plasmid are created we must set the field "timeSpentInHost" of the plasmid
			// as being the time the cell was born (the currentime)
			aPlasmid.timeSpentInHost = SimTimer.getCurrentTime();


		} catch (Exception e) {
			System.out.println("at EpiBac: addPlasmid error " + e.getMessage());
		}
	}

	public boolean sendPlasmid(MultiEpisome aPlasmid, MultiEpiBac partner, double elapsedHGTtime) {

		receivePlasmid(aPlasmid, partner);
		aPlasmid.givePlasmid(elapsedHGTtime);

		return true;
	}

	/**
	 * modified by sonia:
	 * 
	 * In this method we are also recording the information regarding the location and geneology/family/generation
	 * of the recipient cell which will then be written in the agentState.xml output file.
	 * 
	 * @param aPlasmid
	 * @param partner
	 * @return
	 */

	public boolean receivePlasmid(MultiEpisome aPlasmid, MultiEpiBac partner) {
		try {			

			//int plListSize = partner._plasmidHosted.size();
			//ArrayList<String> plHostedNames = new ArrayList<String>();

			// Create a new instance
			MultiEpisome baby = aPlasmid.sendNewAgent();
			baby.setHost(partner);
			//	baby.lastReception = elapsedHGTtime;
			baby.nbhList.clear();
			baby.setDefaultCopyNumber();	
			//baby._newT = SimTimer.getCurrentIter();

			//sonia: 08-06-09
			//update the time this plasmid has entered the recipient, will be used to calculate its fitness cost
			baby.timeSpentInHost = SimTimer.getCurrentTime();

			// register the plasmid to the host	
			//
			partner._plasmidHosted.add(baby);	
			partner.addPlasmidReaction(aPlasmid);

			//StringBuffer partnerInfo = new StringBuffer();
			//partnerInfo.append(partner.getName()+
			//		","+ partner._family + "," + partner._genealogy + 
			//		","+ partner._generation + ","+ partner._location);
			//partnerVector.add(partnerInfo.toString());
			//plasmidVector.add(aPlasmid.getName());
			conjResult = true;
			conjugationEvents++;



		} catch (Exception e) {
			System.out.println("At EpiBac receivePlasmid, exception is..." + e);
			utils.LogFile.writeLog("Error met in EpiBac.receivePlasmid() " +e);		
		}

		if(conjResult){
			System.out.println("ONE TRANSFER!");
		}

		return conjResult;
	}

	public void losePlasmid(MultiEpisome aPlasmid) {
		for (int aReaction : aPlasmid.reactionActive)
			removeReaction(allReactions[aReaction]);
	}

	/**
	 * Add active reaction coded on the plasmid to active reaction of the host
	 * @param aPlasmid
	 */
	public void addPlasmidReaction(MultiEpisome aPlasmid) {
		for (int aReaction : aPlasmid.reactionActive) {
			addActiveReaction(allReactions[aReaction], true);
		}
	}

	public MultiEpiBacParam getSpeciesParam() {
		return (MultiEpiBacParam) _speciesParam;
	}

	/**
	 * Used to write povray files
	 * @return color of the host if empty, color of the first hosted plasmid
	 * else
	 */
	public Color getColor() {

		MultiEpiBacParam param = getSpeciesParam();
		boolean test = false;

		if (_plasmidHosted.size()==0){
			//localAgent.setLocalAgentColor(param.rColor);
			return param.rColor;
		}
		loopCol:
			for(int i=0; i< _plasmidHosted.size(); i++){
				test = _plasmidHosted.get(i).isTransConjugant();
				if(test){
					break loopCol;
				}
			}

		if (test){
			//localAgent.setLocalAgentColor(param.tColor);
			return param.tColor;

		}else{
			//	localAgent.setLocalAgentColor(param.dColor);
			return param.dColor;
		}

		//return localAgent.getLocalAgentColor();
	}



	public double getMuStar() {
		return Math.pow(getNetGrowth()/getSpeciesParam().KSat, getSpeciesParam().nBind);
	}


	public int getPlasmidIndex(String plasmidName, Vector <MultiEpisome> plasmidList){
		int index=0;

		for (int i=0; i<plasmidList.size(); i++){
			if(plasmidList.get(i).getName().equals(plasmidName)){
				index = i;
			}
		}		
		return index;
	}

	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());

		//for(Episome anEpi : _plasmidHosted){	
		tempString.append(",");
		tempString.append("plasmid,copyNumber,transconjugant");

		//}

		return tempString.toString();
	}

	public String writeOutput() {
		//super.writeOutput();	
		//tempString.delete(tempString.length()-2, tempString.length());
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());

		//if (!_plasmidHosted.isEmpty()){
		//tempString.append(",");
		//}


		for(MultiEpisome anEpi : _plasmidHosted){	
			tempString.append(",");
			tempString.append(anEpi.getSpeciesParam().plasmidName + ",");
			tempString.append(anEpi.getCopyNumber()+ ",");
			//sonia:
			//count cells that carry a certain type of plasmid; the copy number of the plasmid is irrelevant
			//for the time being --> that's done in the matlab scripts analysing the agent_sum xml files 

			if(anEpi.isTransConjugant()){
				tempString.append("1");
				//tempString.append( ",");
			}else{
				tempString.append("0");
				//tempString.append( ",");
			}

		}


		//sonia: this will be used when we want to know the information about the partner, to be able to see
		//where the conjugation events take place spatially
		/* if(conjResult){					
			for (int i=0; i< plasmidVector.size(); i++){
				tempString.append("conjugation,plasmid,partnerInfo;");
				tempString.append("1");
				tempString.append("," + plasmidVector.get(i)+"," + partnerVector.get(i));
				tempString.append(";\n");

			}

		}*/
		plasmidVector.clear();
		partnerVector.clear();

		return tempString.toString();
	}

	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		MultiEpiBacParam param = getSpeciesParam();

		fr.write("#declare "+_species.speciesName+"_d = color rgb < ");
		fr.write(((float) param.dColor.getRed()) / 255.0 + " , ");
		fr.write(((float) param.dColor.getGreen()) / 255.0 + " , ");
		fr.write(((float) param.dColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");

		fr.write("#declare "+_species.speciesName+"_r = color rgb < ");
		fr.write(((float) param.rColor.getRed()) / 255.0 + " , ");
		fr.write(((float) param.rColor.getGreen()) / 255.0 + " , ");
		fr.write(((float) param.rColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");

		fr.write("#declare "+_species.speciesName+"_t = color rgb < ");
		fr.write(((float) param.tColor.getRed()) / 255.0 + " , ");
		fr.write(((float) param.tColor.getGreen()) / 255.0 + " , ");
		fr.write(((float) param.tColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}

}