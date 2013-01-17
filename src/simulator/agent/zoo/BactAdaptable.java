/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ___________________________________________________________________________
 * BactAdaptable: a species that can change its reaction based on local conditions 
 * 
 */

/**
 * 
 * @since Nov 2008
 * @version 1.0
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.agent.LocatedAgent;
import simulator.geometry.Bulk;
import simulator.geometry.boundaryConditions.AllBC;

//import utils.LogFile;
import utils.XMLParser;

public class BactAdaptable extends BactEPS {

	// false for off, true for on
	protected Boolean switchState;

	// flag for whether we need to change the state of the switch
	protected Boolean turnSwitchOn = false;
	protected Boolean turnSwitchOff = false;

	// used in conjunction with lag parameters
	protected double timeOfRequestToSwitchOn;
	protected double timeOfRequestToSwitchOff;


	public BactAdaptable() {
		super();
		_speciesParam = new BactAdaptableParam();
	}

	/**
	 * Initialises the progenitor
	 * (This code borrowed from Bacterium class)
	 */
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) {
		// Initialisation of the Bacterium
		super.initFromProtocolFile(aSim, aSpeciesRoot);

		// because the switch parameters are common to a species, they
		// are read in by the Param class
		
		// initially put switch in the 'off' state
		setSwitchToOff();

		// set this so that we can switch reset the switch at the beginning if needed
		timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
		timeOfRequestToSwitchOn = 0;
	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData) {
		// find the position to start at by using length and number of values read
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;

		// read in info from the result file IN THE SAME ORDER AS IT WAS OUTPUT

		// switch state parameters
		//switchState   = Boolean.parseBoolean(singleAgentData[iDataStart]);
		//turnSwitchOn  = Boolean.parseBoolean(singleAgentData[iDataStart+1]);
		//turnSwitchOff = Boolean.parseBoolean(singleAgentData[iDataStart+2]);
		int val;
		val = Integer.parseInt(singleAgentData[iDataStart]);
		if (val==1) switchState = true;
		else		switchState = false;
		val = Integer.parseInt(singleAgentData[iDataStart+1]);
		if (val==1) turnSwitchOn = true;
		else		turnSwitchOn = false;
		val = Integer.parseInt(singleAgentData[iDataStart+2]);
		if (val==1) turnSwitchOff = true;
		else		turnSwitchOff = false;
		timeOfRequestToSwitchOn  = Double.parseDouble(singleAgentData[iDataStart+3]);
		timeOfRequestToSwitchOff = Double.parseDouble(singleAgentData[iDataStart+4]);

		// now go up the hierarchy with the rest of the data
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls
	 */
	protected void internalStep() {
		// check whether we will need to change the switch state
		respondToConditions();
		
		// this is where the switch state is changed once the lag period has passed
		updateActiveReactions();
		
		// once the reactions are set up, everything else goes as normal
		super.internalStep();
	}

	/**
	 * Look in the local region and set whether to change the switch state
	 * (but it will actually be switched later)
	 * 
	 * the switch may be in two states: ON or OFF, and when in each state may be in an
	 * additional two states: waiting to switch, or not waiting to switch. While we are
	 * in the waiting-to-switch state, we treat the switch as if it has already switched
	 * for purposes of testing local conditions (i.e. the planned switch will be cancelled
	 * if the conditions go back to a state not requiring a switch to occur).
	 * 
	 * The four states are:
	 * ON and staying ON
	 * ON and waiting to switch to OFF
	 * OFF and staying OFF
	 * OFF and waiting to switch to ON
	 */
	public void respondToConditions() {

		double localValue=0;

		if(Simulator.isChemostat){
			if (getSpeciesParam().switchType.equals("solute")) {
			
				for (AllBC aBC : _agentGrid.domain.getAllBoundaries()){
					if (aBC.hasBulk()){
						Bulk aBulk = aBC.getBulk();
							if(aBulk.getName().equals("chemostat")){
								localValue = aBulk.getValue(getSpeciesParam()._soluteList[getSpeciesParam().switchControlIndex].soluteIndex);
							}
					}	
				}
				
				
			} else {
				// biomass
				localValue = getParticleMass(getSpeciesParam().switchControlIndex);
			}
			
		}else{
			// get the value needed to check the switch
			if (getSpeciesParam().switchType.equals("solute")) {
				localValue = getSpeciesParam().
						_soluteList[getSpeciesParam().switchControlIndex].
							getValueAround((LocatedAgent) this);
			} else {
				// biomass
				localValue = getParticleMass(getSpeciesParam().switchControlIndex);
			}
		
		}
		
		if ((switchIsOn() && !turnSwitchOff) || (!switchIsOn() && turnSwitchOn)) {
			// state is one of:
			// 		-ON and staying ON
			// 		-OFF and waiting to switch to ON
			// so test whether we need to turn it OFF or cancel the switch
			// (here the test is OPPOSITE of label)
			
			if (getSpeciesParam().switchCondition.equals("lessThan")) {
				if (localValue >= getSpeciesParam().switchValue) {
					if (switchIsOn()) {
						turnSwitchOff = true;
						timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
					} else {
						turnSwitchOn = false;
					}
				}
			} else {
				if (localValue <= getSpeciesParam().switchValue) {
					if (switchIsOn()) {
						turnSwitchOff = true;
						timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
					} else{
						turnSwitchOn = false;
					}
				}
			}
		} else {
			// state is one of:
			// 		-ON and waiting to switch to OFF
			// 		-OFF and staying OFF
			// so test whether we need to turn it ON or cancel the switch
			
			if (getSpeciesParam().switchCondition.equals("lessThan")) {
				if (localValue < getSpeciesParam().switchValue) {
					if (switchIsOn()) {
						turnSwitchOff = false;
					} else {
						turnSwitchOn = true;						
						timeOfRequestToSwitchOn = SimTimer.getCurrentTime();
					}
				}
			} else {
				if (localValue > getSpeciesParam().switchValue) {
					if (switchIsOn()) {
						turnSwitchOff = false;
					} else {
						turnSwitchOn = true;
						timeOfRequestToSwitchOn = SimTimer.getCurrentTime();
					}
				}
			}
		}
	}
	
	public void updateActiveReactions() {
		// leave if we don't need to change the state right now
		if (switchIsOn() && !turnSwitchOff) return;
		if (!switchIsOn() && !turnSwitchOn) return;
		
		// also don't switch if we're still within the lag
		if (!readyToSwitch()) return;


		if (switchIsOn()) {
			setSwitchToOff();
			turnSwitchOff = false;
			
		} else {
			setSwitchToOn();
			turnSwitchOn = false;
			
		}
	}

	// this changes the switch state by turning reactions on or off
	public void setSwitchToOn() {
		// turn off the reactions that were previously on
		for (int aReac : getSpeciesParam().offStateReactions) {
			switchOffreaction(allReactions[aReac]);
		}

		// turn on the reactions that should now be on
		for (int aReac : getSpeciesParam().onStateReactions) {
//			
			switchOnReaction(allReactions[aReac]);
		}

		setSwitchState(true);
		
		
	}
	
	// this changes the switch state by turning reactions on or off
	public void setSwitchToOff() {

		// turn off the reactions that were previously on
		for (int aReac : getSpeciesParam().onStateReactions) {
			switchOffreaction(allReactions[aReac]);
		}

		// turn on the reactions that should now be on
		for (int aReac : getSpeciesParam().offStateReactions) {
			switchOnReaction(allReactions[aReac]);
		}

		setSwitchState(false);
	
	}

	// returns false when off, true when on
	public Boolean switchIsOn() {
		return switchState;
	}
	
	// set the state to true or false
	public void setSwitchState(Boolean newState) {
		switchState = newState;
	}
	
	// return true or false if ready to switch
	// this prevents switching during the initial adaptation lag period
	public Boolean readyToSwitch() {
		// don't switch if we're still within the lag

		// TODO: make this a distribution around the time rather than discrete decision
		if (switchIsOn()) {
			// if it's on already, need to consider off-state lag
			if (SimTimer.getCurrentTime() - timeOfRequestToSwitchOff < getSpeciesParam().lagSwitchOff)
				return false;
		} else {
			// if it's off, need to consider on-state lag
			if (SimTimer.getCurrentTime() - timeOfRequestToSwitchOn < getSpeciesParam().lagSwitchOn)
				return false;
		}

		return true;
	}

	public BactAdaptableParam getSpeciesParam() {
		return (BactAdaptableParam) _speciesParam;
	}
	

	/* _______________ FILE OUTPUT _____________________ */

	public String sendHeader() {
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",");
		
		// switch state and timing info
		tempString.append("state,turnOn,turnOff,timeRequestTurnOn,timeRequestTurnOff");

		return tempString.toString();
	}

	public String writeOutput() {
		// write the data matching the header file
		StringBuffer tempString = new StringBuffer(super.writeOutput());
		tempString.append(",");

		// switch state and timing info

		if (switchState)   tempString.append("1,");
		else			   tempString.append("0,");
		if (turnSwitchOn)  tempString.append("1,");
		else			   tempString.append("0,");
		if (turnSwitchOff) tempString.append("1,");
		else			   tempString.append("0,");
		tempString.append(timeOfRequestToSwitchOn+","+timeOfRequestToSwitchOff);

		return tempString.toString();
	}
	
	
	/**
	 * Used to write povray files (replaces the version in LocatedAgent)
	 */
	public String getName() {
		if (switchIsOn())
			return _species.speciesName+"_ON";

		return _species.speciesName+"_OFF";
	}
/**
	 * Used to write povray files
	 */
	public Color getColor() {
		if (switchIsOn())
			return getSpeciesParam().onColor;

		return getSpeciesParam().offColor;
	}

	/**
	 * this writes a color definition to the passed-in file; meant for later use in macros
	 * This overrules the version in SpecialisedAgent
	 * 
	 * @param theFile
	 */
	public void writePOVColorDefinition(FileWriter fr) throws IOException {
		BactAdaptableParam param = getSpeciesParam();
		
		fr.write("#declare "+_species.speciesName+"_ON = color rgb < ");
		fr.write(((float) param.onColor.getRed()) / 255.0 + " , ");
		fr.write(((float) param.onColor.getGreen()) / 255.0 + " , ");
		fr.write(((float) param.onColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
		
		fr.write("#declare "+_species.speciesName+"_OFF = color rgb < ");
		fr.write(((float) param.offColor.getRed()) / 255.0 + " , ");
		fr.write(((float) param.offColor.getGreen()) / 255.0 + " , ");
		fr.write(((float) param.offColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}
}
