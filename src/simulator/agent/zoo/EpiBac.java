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

import idyno.SimTimer;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import org.jdom.Element;
import utils.ExtraMath;
import utils.XMLParser;
import utils.LogFile;
import simulator.Simulator;
import simulator.agent.*;
import simulator.geometry.ContinuousVector;

// bvm 30.1.2009
// made this derive from BactEPS, not Bacterium
public class EpiBac extends BactEPS {

	/* Parameters mutated from species parameters ___________________________ */

	/* Parameters specific to the agent _____________________________________ */
	// Plasmid hosted by this agent
	private LinkedList<Episome> _plasmidHosted = new LinkedList<Episome>();
	private double _lastReception=0;
	private double _lastExchange=0;
	private int _nCopy=0;
	private int _status=0;

	/* _________________________ CONSTRUCTOR _____________________________ */
	/**
	 * Empty constructor ; called to build a progenitor ; the speciesParameter
	 * have to be defined later
	 */
	public EpiBac() {
		super();
		_speciesParam = new EpiBacParam();
	}

	public Object clone() throws CloneNotSupportedException {
		EpiBac o = (EpiBac) super.clone();
		o._plasmidHosted = new LinkedList<Episome>();
		Episome newEpisome;
		for (Episome anEpisome : _plasmidHosted) {
			newEpisome = (Episome) anEpisome.clone();
			newEpisome.setHost(o);
			o._plasmidHosted.add(newEpisome);
		}
		return (Object) o;
	}

	/**
	 * Called during species creation to build the progenitor
	 */
	public void initFromProtocolFile(Simulator aSimulator,
			XMLParser aSpeciesRoot) {
		// Initialisation of the Located agent
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);

		// Create hosted plasmids
		for (Element aSpeciesMarkUp : aSpeciesRoot.buildSetMarkUp("plasmid"))
			addPlasmid(aSpeciesMarkUp.getAttributeValue("name"));

		// Genealogy and size management
		init();
	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// find the position to start at by using length and number of values read
		int nValsRead = 4;
		int iDataStart = singleAgentData.length - nValsRead;

		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT

		// HGT parameters
		_status        = Integer.parseInt(singleAgentData[iDataStart]);
		_nCopy         = Integer.parseInt(singleAgentData[iDataStart+1]);
		_lastReception = Double.parseDouble(singleAgentData[iDataStart+2]);
		_lastExchange  = Double.parseDouble(singleAgentData[iDataStart+3]);

		// now go up the hierarchy with the rest of the data
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	//	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
	//		super.initFromResultFile(aSim,singleAgentData);
	//
	//		// to make the length-i approach work, the singleAgentData string MUST have the 
	//		// relevant values at the end
	//		// this means any derived class must strip its own values from the array before calling
	//		// this routine. For now, just use the hard-coded position
	//		
	//		//int len = singleAgentData.length;
	//		//_status = Integer.parseInt(singleAgentData[len-4]);
	//		//_nCopy = Integer.parseInt(singleAgentData[len-3]);
	//		//_lastReception = Double.parseDouble(singleAgentData[len-2]);
	//		//_lastExchange = Double.parseDouble(singleAgentData[len-1]);
	//		
	//		// family,genealogy,generation,locationX,locationY,locationZ,
	//		// biomass,inert,capsule,
	//		// radius,totalRadius,
	//		// growthRate,volumeRate,birthday,
	//		// status,copyNumber,lastReception,lastExchange
	//		_status = Integer.parseInt(singleAgentData[14]);
	//		_nCopy = Integer.parseInt(singleAgentData[15]);
	//		_lastReception = Double.parseDouble(singleAgentData[16]);
	//		_lastExchange = Double.parseDouble(singleAgentData[17]);
	//	}

	public EpiBac sendNewAgent() throws CloneNotSupportedException {
		EpiBac baby = (EpiBac) this.clone();
		baby.init();
		return baby;
	}

	public void createNewAgent(ContinuousVector position) {
		super.createNewAgent(position);
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
		EpiBac baby = (EpiBac) sendNewAgent();
		baby.mutateAgent();

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
		for (int i = 0; i < _plasmidHosted.size(); i++) {
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
		updateSize();

		// test if the EPS capsule has to be excreted
		manageEPS();

		// Test division and death treshold
		if (willDivide()) divide();
		if (willDie()) die(true);		

		// Try to conjugate
		try{
			conjugate();
			updateStatus();
		}catch(Exception e){
			LogFile.writeError("Failed to conjugate", "EpiBac.internalStep()");
		}
	}

	/**
	 * Remove a plasmid if its copy number reaches 0
	 */
	protected void checkMissingPlasmid() {
		int nIter = _plasmidHosted.size();
		Episome anEpi;
		for(int i=0;i<nIter;i++){
			anEpi = _plasmidHosted.removeFirst();
			if (anEpi.getCopyNumber() <= 0)
				losePlasmid(anEpi);
			else
				_plasmidHosted.addLast(anEpi);				
		}
	}

	/**
	 * test if this cell can be a recipient for a given plasmid
	 * @param aPlasmid
	 * @return true if this cell is compatible
	 */
	public boolean isCompatible(Episome aPlasmid) {
		boolean out = true;
		for (Episome hostPlasmid : _plasmidHosted)
			out &= hostPlasmid.isCompatible(aPlasmid);

		return out;
	}

	/**
	 * Search a recipient in your nbh, and try to initiate a conjugation with
	 * him TODO : cells are attracted
	 */
	public void conjugate() {
		// For each plasmid ready to conjugate search a partner and conjugate
		for (Episome aPlasmid : _plasmidHosted)
			if (aPlasmid.isReadyToConjugate())
				searchConjugation(aPlasmid);
	}

	/**
	 * Remove agent and all references from the system
	 */
	public void die(boolean isStarving) {
		super.die(isStarving);
	}

	/* __________________ CONJUGATION ___________________________ */

	/**
	 * Search a partner and try to send him a plasmid
	 * 
	 * @param aPlasmid
	 * @return issue of the transfer
	 */
	public boolean searchConjugation(Episome aPlasmid) {
		LocatedAgent aLoc;
		Boolean hasTransfered = false;

		// Build a neighborhood including only Bacteria
		buildNbh(aPlasmid.getPilusRange());
		if(_myNeighbors.isEmpty()) return false;

		// Find a partner and try to send him a plasmid
		int iLoop = 0;
		double maxTest = getMaxTest();
		while (iLoop < maxTest) {			
			aLoc = pickNeighbor();			
			iLoop++;
			if (aLoc instanceof EpiBac)
				hasTransfered = tryToSendPlasmid(aPlasmid, (EpiBac) aLoc);
		}
		return hasTransfered;
	}

	/**
	 * List all cells in a given nbh : at the end of the method, the field
	 * listNbh contains all clocatedAgents located in the nbh
	 * 
	 * @param nbhRadius
	 */
	public void buildNbh(double nbhRadius) {
		// Manhattan perimeter
		getPotentialShovers(nbhRadius);

		// Now remove too far agents (apply circular perimeter)
		LocatedAgent aLocAgent;
		for (int iter = 0; iter < _myNeighbors.size(); iter++) {
			aLocAgent = _myNeighbors.removeFirst();
			if (aLocAgent == this)
				continue;
			if (!(aLocAgent instanceof Bacterium))
				continue;

			// getDistance(aLocAgent) gets distance between cell centers
			// nbhRadius is the pilus length, and should give the distance OUTSIDE
			// the donor agent that touches a recipient agent, and so we need to 
			// subtract the radii from getDistance() or add the radii to nbhRadius 
			if (getDistance(aLocAgent) < (nbhRadius + aLocAgent.getRadius(false) + this.getRadius(false)) )
				_myNeighbors.addLast(aLocAgent);
		}
		Collections.shuffle(_myNeighbors, ExtraMath.random);
	}

	/**
	 * Test transfer proficiency and if success, send the plasmid
	 * 
	 * @param aPlasmid
	 * @param aTarget
	 * @return issue of the test
	 */
	public boolean tryToSendPlasmid(Episome aPlasmid, EpiBac aTarget) {
		if (aTarget.isCompatible(aPlasmid) && aPlasmid.testProficiency()) {
			//LogFile.writeLog("transfer "+aPlasmid.getHost().sendName()+"->"+aTarget.sendName());

			// 19.3.09: bvm add more descriptive text to log
			EpiBac aHost = aPlasmid.getHost();
			String transType = "Unknown transfer ";

			// _status == 0: recipient
			// _status == 1: donor
			// _status == 11: transconjugant able to transfer
			// _status == 101: donor has transferred
			// _status == 111: transconjugant has transferred
			int hostStatus = aHost.getStatus();
			if (hostStatus == 1) {
				transType = "Donor initial transfer ";
			} else if (hostStatus == 101) {
				transType = "Donor subsequent transfer ";
			} else if (hostStatus == 11) {
				transType = "Transconjugant initial transfer ";
			} else if (hostStatus == 111) {
				transType = "Transconjugant subsequent transfer ";
			}

			LogFile.writeLog(transType+aHost.sendName()+" -> "+aTarget.sendName()+
					" at time "+SimTimer.getCurrentTime());
			return aTarget.receivePlasmid(aPlasmid);
		} else
			return false;
	}

	/* _______________________ HIGH LEVEL METHOD ____________________________ */
	/**
	 * Add a new plasmid to the list of hosted plasmids ; based on the
	 * speciesname of the plasmid
	 */
	public boolean addPlasmid(String plasmidName) {
		try {
			Episome aPlasmid = (Episome) _species.getSpecies(plasmidName)
			.sendNewAgent();

			// this should match the receivePlasmid() routine
			aPlasmid.setHost(this);
			_plasmidHosted.add(aPlasmid);
			addPlasmidReaction(aPlasmid);
			aPlasmid.givePlasmid(aPlasmid);

			return true;
		} catch (Exception e) {
			utils.LogFile.writeLog("Error met in EpiBac.addPlasmid()");
			return false;
		}
	}

	public boolean receivePlasmid(Episome aPlasmid) {
		try {
			// Create a new instance
			Episome baby = aPlasmid.sendNewAgent();

			// Set new plasmid descriptors (copy number, host, timers)
			baby.setHost(this);

			// register the plasmid to the host
			_plasmidHosted.add(baby);
			addPlasmidReaction(baby);
			aPlasmid.givePlasmid(baby);

			return true;
		} catch (Exception e) {
			utils.LogFile.writeLog("Error met in EpiBac.receivePlasmid()");
			return false;
		}
	}

	public void losePlasmid(Episome aPlasmid) {		
		for (int aReaction : aPlasmid.reactionActive)
			removeReaction(allReactions[aReaction]);
		aPlasmid.die();
		LogFile.writeLog("Plasmid lost "+sendName());
	}

	/**
	 * Add active reaction coded on the plasmid to active reaction of the host
	 * 
	 * @param aPlasmid
	 */
	public void addPlasmidReaction(Episome aPlasmid) {
		for (int aReaction : aPlasmid.reactionActive) {
			addActiveReaction(allReactions[aReaction], true);
		}
	}

	public EpiBacParam getSpeciesParam() {
		return (EpiBacParam) _speciesParam;
	}

	public double getMaxTest() {
		// ======== no growth dependence ======== 
		//		return ((EpiBacParam) _speciesParam).scanSpeed*SimTimer.getCurrentTimeStep();


		// ======== all in one ========
		// return 0 when below lowTonus
		// return 1 when above highTonus
		// interpolate in-between
		//
		// (for no growth dependence, set highTonus = -Double.MAX_VALUE)
		// (for step dependence, set high & low tonus the same)

		double lowTonus = ((EpiBacParam) _speciesParam).lowTonusCutoff;
		double highTonus = ((EpiBacParam) _speciesParam).highTonusCutoff;
		double theTonus = sendTonus();

		if (theTonus >= highTonus) {
			//			System.out.println("high case: "+theTonus);
			// high tonus, so return maximum (same effect as no growth dependence)
			return ((EpiBacParam) _speciesParam).scanSpeed*
			SimTimer.getCurrentTimeStep();
		} else if (theTonus < lowTonus) {
			//			System.out.println("low case: "+theTonus);
			// too low, so return 0
			return 0.;
		} else {
			//			System.out.println("middle case: "+theTonus);
			// middle case, so do linear interpolation
			double vs = ((EpiBacParam) _speciesParam).scanSpeed;
			vs = vs*(theTonus-lowTonus)/(highTonus-lowTonus);
			return vs*SimTimer.getCurrentTimeStep();
		}
	}


	/**
	 * Used to write povray files (replaces the version in LocatedAgent)
	 */
	public String getName() {

		if (_plasmidHosted.size() == 0)
			return _species.speciesName+"_r";

		boolean isTransHost = _plasmidHosted.getFirst().lastReception > _birthday;
		if (isTransHost)
			return _species.speciesName+"_t";
		else
			return _species.speciesName+"_d";
	}

	/**
	 * Used to write povray files
	 */
	public Color getColor() {
		EpiBacParam param = getSpeciesParam();

		// recipients have no plasmid
		if (_plasmidHosted.size() == 0)
			return param.rColor;

		// transconjugant received the plasmid after birth;
		// donor received the plasmid before/at birth
		boolean isTransHost = _plasmidHosted.getFirst().lastReception > _birthday;
		if (isTransHost)
			return param.tColor;
		else
			return param.dColor;
	}

	/* _______________ FILE OUTPUT _____________________ */

	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",");

		tempString.append("status,copyNumber,lastReception,lastExchange");

		return tempString.toString();
	}

	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());
		tempString.append(",");

		// what about _plasmidHosted?
		tempString.append(_status+","+_nCopy+","+_lastReception+","+_lastExchange);

		return tempString.toString();
	}

	//	public String sendHeader() {
	//		super.sendHeader();
	//		// tempString.delete(tempString.length()-2, tempString.length());
	//		tempString.append(",status,copyNumber,lastReception,lastExchange");
	//		return tempString.toString();
	//	}
	//
	//	public String writeOutput() {
	//		super.writeOutput();
	//
	//		// remove the semicolon and endline character that were present
	//		// so that we can append the extra information
	//		tempString.delete(tempString.length() - 2, tempString.length());
	//
	//		tempString.append(",");
	//		tempString.append(_status);
	//		tempString.append(",");
	//		tempString.append(_nCopy);
	//		tempString.append(",");
	//		tempString.append(_lastReception);
	//		tempString.append(",");
	//		tempString.append(_lastExchange);
	//		tempString.append(";\n");
	//		return tempString.toString();
	//	}

	public double sendTonus(){
		// this returns the net growth rate as a fraction of the maximum rate,
		// and the value is as large as 1 but may also be negative

		double tonusMax=0;
		int reacIndex;

		for (int iReac = 0; iReac<reactionActive.size(); iReac++) {
			// Compute the growth rate
			reacIndex = reactionActive.get(iReac);
			tonusMax += reactionKinetic[reacIndex][0]*particleYield[reacIndex][0];			
		}

		tonusMax *= particleMass[0];
		return _netGrowthRate/tonusMax;
	}

	public void updateStatus(){
		_nCopy=0;
		_status=0;

		// _status == 0: recipient
		// _status == 1: donor
		// _status == 11: transconjugant able to transfer
		// _status == 101: donor has transferred
		// _status == 111: transconjugant has transferred


		// ones digit specifies number of plasmids hosted
		for(Episome anEpi:_plasmidHosted){
			_lastReception = Math.max(_lastReception,anEpi.lastReception);
			_lastExchange = Math.max(_lastExchange,anEpi.lastExchange);
			_nCopy = Math.max(_nCopy,anEpi.getCopyNumber());
			_status+=1;
		}
		if(_nCopy==0)
			_status=0;
		else {
			// transconjugant received the plasmid after birth;
			// donor received the plasmid before/at birth

			// tens digit specifies if it is a transconjugant
			if(_lastReception > _birthday) _status += 10;

			// hundreds digit specifies whether it has transferred a plasmid
			if(_lastExchange > _birthday) _status += 100;
		}
	}

	// return the current status for use elsewhere
	// _status == 0: recipient
	// _status == 1: donor
	// _status == 11: transconjugant able to transfer
	// _status == 101: donor has transferred
	// _status == 111: transconjugant has transferred
	public int getStatus() {
		return _status;
	}

	/**
	 * this writes a color definition to the passed-in file; meant for later use in macros
	 * This overrules the version in SpecialisedAgent
	 * 
	 * @param theFile
	 */
	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		EpiBacParam param = getSpeciesParam();

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