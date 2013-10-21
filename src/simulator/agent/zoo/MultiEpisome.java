/**
 * Project MicoDyna
 * Class describing the dynamic of one kind of plasmid in a cell
 * 
 * @since June 2006
 * @version 0.1
 * @author Laurent Lardon, DTU Denmark
 */

package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.LinkedList;

import org.jdom.Element;

import utils.XMLParser;
import utils.ExtraMath;

import idyno.SimTimer;
import simulator.agent.*;
import simulator.reaction.Reaction;
import simulator.Simulator;

public class MultiEpisome extends InfoAgent {

	// Serial version used for the serialisation of the class
	private static final long   serialVersionUID  = 1L;

	//sonia:I've changed it to public
	public MultiEpiBac              _host;

	/* Division management ___________________________________ */
	
	//sonia:I've changed it to public
	public double              _nCopy;

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

		return (Object) o;
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
		
		for (Element aReactionMarkUp : xmlMarkUp.buildSetMarkUp("reaction")) {
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));

			// Add the reaction to the list of known (and active) reactions
			reactionKnown.add(reacIndex);
			if (aReactionMarkUp.getAttributeValue("status").equals("active")) {
				reactionActive.add(reacIndex);
			}
		}
		
	}

	public MultiEpisome reset(MultiEpisome aPlasmid){
		aPlasmid._generation = 0;
		aPlasmid._genealogy = 0;
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
		_genealogy = 0;
		//sonia 01.03.2010 changed from -1 to 0
		lastExchange = 0;
		lastReception = SimTimer.getCurrentTime();
	}

	/**
     * 
     */
	public MultiEpisome sendNewAgent() throws CloneNotSupportedException {
		MultiEpisome baby = (MultiEpisome) this.clone();
		baby.init();
		return baby;
	}

	/**
     * 
     */
	public void createNewAgent() {
		try {
			// Clone the plasmid
			MultiEpisome baby = this.sendNewAgent();

			// Mutate parameters
			baby.mutatePop();

			// Register the plasmid (species population)
			baby.registerBirth();

		} catch (CloneNotSupportedException e) {
			System.out.println("at Episome: createNewAgent error " + e);
		}
	}

	public void mutatePop() {
		// Mutate inherited parameters
		super.mutatePop();

		// Now mutate your parameters
	}

	public void registerBirth() {
		_species.notifyBirth();
	}

	/* __________________________ CELL DIVISION ____________________________ */

	public void makeKid() throws CloneNotSupportedException {
		// Clone the plamid
		MultiEpisome baby = this.sendNewAgent();

		// Mutate parameters
		baby.mutateAgent();

		// Register the plasmid (species population)
		baby.registerBirth();
	}

	public void mutateAgent() {
		// Mutate inherited parameters
		super.mutateAgent();
		// lastReception = -(getSpeciesParam().receptionLag +1);
		// lastExchange = -(getSpeciesParam().exchangeLag +1);
		//sonia 11.10.2010
		lastReception=0;
		lastExchange=0;
		// Now mutate your parameters
	}

	/* _______________________________________________________________________ */

	public void internalStep() {
		//sonia 6.10.09
		// for now we will only have a fixed number of plasmid copies, hence there will be 
		//no plasmid replication
		//replication();
	}


	/**
     * Test if the episome can be transfered
     */
	public boolean isReadyToConjugate(double elapsedHGTtime) {
		
		//EpisomeParam param = getSpeciesParam();
		//isHot = false;
		boolean conjugate = false;
		//sonia 5.10.2010
	//	boolean exchange = false;
		//boolean receive = false;
		
		if (_nCopy < 1){

			}else{
			conjugate = true;
		
		}
	
		return conjugate;
/*		
		if(Simulator.isChemostat){
			//in the chemostat situation the plasmid is always ready to conjugate, for now.
			return true;
		
		}else{
		
		//sonia 4.10.2010
			
		//	System.out.println("checking if plasmid is ready to conjugate");
			
		// You have given a plasmid a few minutes ago
		if (elapsedHGTtime-this.lastExchange>param.exchangeLag){
			exchange = true;
		}

		if ((elapsedHGTtime-this.lastReception>param.receptionLag)){
			receive =true;
		}
		
		     if (exchange & receive){
		    //	 System.out.println("true");
		     	  return true;
		     }else{
		    //	 System.out.println("false");
		    	  return false;
		     }
		}
		*/
		//return receive;

	}



	/**
	 * NOT BEING USED
     * Number of plasmid copies is doubled at the doubling speed rate until the
     * maximal value is reached
     */
	public void replication() {
		int nCopyMax = getSpeciesParam().nCopyMax;

		if (_nCopy>=nCopyMax) {
			_nCopy = nCopyMax;
			return;
		}
		double doublingTime = getSpeciesParam().replicationSpeed;
		double muStar = _host.getMuStar();
		_nCopy += doublingTime*muStar/(1+muStar)*SimTimer.getCurrentTimeStep()*_nCopy;
	}

	public void die() {
		_species.notifyDeath();
	}

	public double getPilusLength() {
		return _pilusLength;
	}

	/**
     * You are doing a conjugation ! Update your parameters
     */
	public void givePlasmid(double elapsedHGTtime) {
		
		//sonia 4.10.2010
		// Update your time counter for conjugation-able recovery
		lastExchange = elapsedHGTtime;
	}
	
	/**
	 * sonia:14-05-09
	 * Set the lastReception field to the current time which corresponds to the time
	 * the cell received the plasmid.
	 */
	public void receivedPlasmid(){
		
		MultiEpisomeParam param = this.getSpeciesParam();
		//lastReception += param.receptionLag;
		//lastReception = SimTimer.getCurrentTime();
	}
	

/*	*//**
     * Called during the division of an EpiBac ; apply the segregation of
     * plasmids, modify the number for the plasmid calling the method and sends
     * the number for the other one
     * @return
     *//*
	public void segregation(Episome aPlasmid) {
		int plCopy=0;
		double nSplit = 0;
		if (getSpeciesParam().isHighCopy) {
			// Simulates a binomial segregation
			for (int iCopy = 0; iCopy<aPlasmid._nCopy; iCopy++) {
				nSplit += (ExtraMath.getNormRand()+2)/2;
			}
			plCopy = ((int) nSplit);
		} else {
			if (ExtraMath.getUniRand()< getSpeciesParam().lossProbability) 
				//plCopy = (int) Math.round(aPlasmid._nCopy/2);
			
				plCopy = (int)aPlasmid._nCopy;
			
			else
				plCopy = 0;
		}
		
		aPlasmid._nCopy = aPlasmid._nCopy - plCopy;
	}*/

	
	/**
	 * Called during the division of an EpiBac ; apply the segregation of
	 * plasmids, modify the number for the plasmid calling the method and sends
	 * the number for the other one
	 * 
	 * @return
	 */
	public void segregation(MultiEpisome aPlasmid) {
		if (ExtraMath.getUniRandDbl() > getSpeciesParam().lossProbability) {
			_nCopy = 1;
			aPlasmid._nCopy = 1;
		} else {
			_nCopy = 0;
			aPlasmid._nCopy = 1;
		}
	}
	
	
	
	
	public MultiEpisomeParam getSpeciesParam() {
		return (MultiEpisomeParam) _speciesParam;
	}

	public int getCopyNumber() {
		return (int) _nCopy;
	}

	public void setCopyNumber(int value) {
		_nCopy = value;
	}
	
/**
 * sonia:
 * to be used in increasing the copy number of a plasmid when this is received by a cell
 * that already contains a plasmid of this type. - which will never happen due to entry exclusion systems...
 */
	public void increaseCopyNumber(){
		_nCopy ++;
	}
	
	
	public void setDefaultCopyNumber() {
		_nCopy = getSpeciesParam().nCopy;
	}

	/**
	 * sonia:
	 * used in conjugation events: when a cell receives a plasmid, it will start off with just one copy
	 */
	public void setInitialCopyNumber(){
		_nCopy = 1;
	}
	
	//sonia: 23-07-09
	// now we can distinguish between vertical and horizontal transmission of the plasmid
	//if the time of its reception is higher than the time at which its host was born, then the plasmid 
	//is harboured by a transconjugant cell
	public boolean isTransConjugant() {
		boolean transconjugant = false;
		if (this.lastReception > this._host.getBirthday()){
			transconjugant = true;
		}
		return transconjugant;
		//return SimTimer.getCurrentTime()-lastReception<getSpeciesParam().receptionLag;
	}
	



	// TODO
	public String sendHeader() {
		return "locationX;locationY;locationZ;mass;radius;growthRate";
	}

	// TODO
	public String writeOutput() {
		// Now send your data
		//tempString.append("copyNumber,transconjugant;");
		return (_species.speciesName);
	}



	//@Override
	//protected void conjugate(double elapsedHGTtime) {
		// TODO Auto-generated method stub
		
	//}
	
	
	
}
