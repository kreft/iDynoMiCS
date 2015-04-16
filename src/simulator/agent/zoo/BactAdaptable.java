/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to
 * store parameters for these agent types.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.agent.zoo;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SoluteGrid;
import utils.XMLParser;

/**
 * \brief Creates a Bacterium agent object that can change its reaction based
 * on local conditions.
 * 
 * The BactAdaptable species derives from the BactEPS species, but adds the
 * ability to change its set of active reactions based on definable
 * conditions. The mechanism for this switching behaviour is an ON/OFF state 
 * switch, with each state having a particular set of active reactions that
 * are set as inactive when the switch is in the opposite state. You may use 
 * local solute concentrations or agent biomass amounts as switch conditions; 
 * this example mark-up uses the local oxygen (MyO2) concentration as the
 * switch condition.
 * 
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of
 * Engineering Sciences and Applied Mathematics, Northwestern University (USA).
 */
public class BactAdaptable extends BactEPS 
{
	/**
	 * Boolean that notes whether this switch is set to be on (true) or false
	 * (off).
	 */
	protected Boolean switchState;
	
	/**
	 * Boolean flag noting whether we need to change the state of the switch
	 * to on.
	 */
	protected Boolean turnSwitchOn = false;
	
	/**
	 * Boolean flag noting whether we need to change the state of the switch
	 * to off.
	 */
	protected Boolean turnSwitchOff = false;
	
	/**
	 * Used in conjunction with lag parameters, notes the time that the
	 * simulation requested the switch be set to on.
	 */
	protected Double timeOfRequestToSwitchOn;
	
	/**
	 * Used in conjunction with lag parameters, notes the time that the
	 * simulation requested the switch be set to off.
	 */
	protected Double timeOfRequestToSwitchOff;
	
	/**
	 * \brief Constructor used to generate progenitor and initialise an object
	 * to store relevant parameters.
	 */
	public BactAdaptable()
	{
		super();
		_speciesParam = new BactAdaptableParam();
	}
	
	/**
	 * \brief Creates a BactAdaptable agent from the parameters specified in
	 * the XML protocol file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aSpeciesRoot	A species mark-up within the specified protocol
	 * file.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) 
	{
		/*
		 * Initialization of the Bacterium.
		 * 
		 * Because the switch parameters are common to a species, they are
		 * read in by the Param class.
		 */
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		/*
		 * Initially put switch in the 'off' state.
		 */
		setSwitchToOff();
		/*
		 * Set this so that we can switch reset the switch at the beginning
		 * if needed.
		 */
		timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
		timeOfRequestToSwitchOn = 0.0;
	}

	/**
	 * \brief Create an agent using information in a previous state or
	 * Initialization file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialization file
	 * that is used to recreate this agent.
	 */
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		/*
		 * Find the position to start at by using length and number of values
		 * read.
		 */
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * Read in info from the result file IN THE SAME ORDER AS IT WAS
		 * OUTPUT.
		 * 
		 * Switch state parameters.
		 */
		switchState = singleAgentData[iDataStart].equals("1");
		turnSwitchOn = singleAgentData[iDataStart+1].equals("1");
		turnSwitchOff = singleAgentData[iDataStart+2].equals("1");
		timeOfRequestToSwitchOn  =
							Double.parseDouble(singleAgentData[iDataStart+3]);
		timeOfRequestToSwitchOff =
							Double.parseDouble(singleAgentData[iDataStart+4]);
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	/**
	 * \brief Called at each time step of the simulation to compute mass
	 * growth and update radius, mass, and volume. In this case also checks
	 * switch state.
	 * 
	 * Called at each time step of the simulation (under the control of the
	 * method Step of the class Agent) to compute mass growth and update
	 * radius, mass, and volume. Also determines whether the agent has reached
	 * the size at which it must divide, and monitors agent death. In this
	 * case also checks switch state.
	 */
	@Override
	protected void internalStep() 
	{
		/*
		 * Check whether we will need to change the switch state.
		 */
		respondToConditions();
		/*
		 * This is where the switch state is changed once the lag period has
		 * passed.
		 */
		updateActiveReactions();
		/*
		 * Once the reactions are set up, everything else goes as normal.
		 */
		super.internalStep();
	}
	
	/**
	 * \brief Examines the local region and determines whether to set the flag
	 * noting the switch should be changed.
	 * 
	 * Look in the local region and set whether to change the switch state
	 * (but it will actually be switched later). The switch may be in two
	 * states: ON or OFF, and when in each state may be in an additional two
	 * states: waiting to switch, or not waiting to switch. 
	 * While we are in the waiting-to-switch state, we treat the switch as if
	 * it has already switched for purposes of testing local conditions (i.e.
	 * the planned switch will be cancelled if the conditions go back to a
	 * state not requiring a switch to occur).
	 * 
	 * The four states are:
	 * A. ON and staying ON;
	 * B. ON and waiting to switch to OFF;
	 * C. OFF and staying OFF;
	 * D. OFF and waiting to switch to ON
	 */
	public void respondToConditions()
	{
		Double localValue = 0.0;
		int sc = getSpeciesParam().switchControlIndex;
		Double switchValue = getSpeciesParam().switchValue;
		Boolean switchLessThan = 
						getSpeciesParam().switchCondition.equals("lessThan");
		/*
		 * If the switch type is solute, find the local concentration.
		 * Otherwise, find the relevant biomass.
		 */
		if ( getSpeciesParam().switchType.equals("solute"))
		{
			SoluteGrid sg = getSpeciesParam()._soluteList[sc];
			if ( Simulator.isChemostat)
			{
				localValue = _agentGrid.domain.getChemostat().
													getValue(sg.soluteIndex);
			}
			else
				localValue = sg.getValueAround(this);
		}
		else
			localValue = getParticleMass(sc);
		/*
		 * Now we've found the local value, use it to determine the response.
		 */
		if ((switchState && !turnSwitchOff) || (!switchState && turnSwitchOn))
		{
			/*
			 * State is one of:
			 * 		-ON and staying ON
			 * 		-OFF and waiting to switch to ON
			 * so test whether we need to turn it OFF or cancel the switch
			 * (here the test is OPPOSITE of label)
			 */
			if ( switchLessThan )
			{
				if ( localValue >= switchValue )
				{
					if ( switchState )
					{
						turnSwitchOff = true;
						timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
					}
					else
						turnSwitchOn = false;
				}
			}
			else
			{
				if ( localValue <= switchValue )
				{
					if ( switchState )
					{
						turnSwitchOff = true;
						timeOfRequestToSwitchOff = SimTimer.getCurrentTime();
					}
					else
						turnSwitchOn = false;
				}
			}
		}
		else
		{
			/*
			 * State is one of:
			 * 		-ON and waiting to switch to OFF
			 * 		-OFF and staying OFF
			 * so test whether we need to turn it ON or cancel the switch
			 */
			if ( switchLessThan )
			{
				if ( localValue < switchValue )
				{
					if ( switchState )
						turnSwitchOff = false;
					else
					{
						turnSwitchOn = true;						
						timeOfRequestToSwitchOn = SimTimer.getCurrentTime();
					}
				}
			}
			else
			{
				if ( localValue > switchValue )
				{
					if ( switchState )
						turnSwitchOff = false;
					else
					{
						turnSwitchOn = true;
						timeOfRequestToSwitchOn = SimTimer.getCurrentTime();
					}
				}
			}
		}
	}
	
	/**
	 * \brief Used by the internal step to update the active reactions
	 * dependent on the state of the switch.
	 * 
	 * Checks to ensure adherence to lag time.
	 */
	public void updateActiveReactions() 
	{
		/*
		 * Leave if we don't need to change the state right now.
		 */
		if ( switchState && ! turnSwitchOff )
			return;
		if ( ! switchState && ! turnSwitchOn )
			return;
		/*
		 * Also don't switch if we're still within the lag.
		 */
		if ( ! readyToSwitch() )
			return;
		/*
		 * 
		 */
		if ( switchState )
		{
			setSwitchToOff();
			turnSwitchOff = false;
		}
		else
		{
			setSwitchToOn();
			turnSwitchOn = false;
		}
	}

	/**
	 * \brief Changes the switch state to on by turning reactions on or off,
	 * dependent on the state the reaction should be in when the switch is on.
	 */
	public void setSwitchToOn()
	{
		/*
		 * Turn off the reactions that were previously on.
		 */
		for ( int aReac : getSpeciesParam().offStateReactions )
			switchOffreaction(allReactions[aReac]);
		/*
		 * Turn on the reactions that should now be on.
		 */
		for ( int aReac : getSpeciesParam().onStateReactions )
			switchOnReaction(allReactions[aReac]);
		/*
		 * 
		 */
		setSwitchState(true);
	}
	
	/**
	 * \brief Changes the switch state to off by turning reactions on or off,
	 * dependent on the state the reaction should be in when the switch is off.
	 */
	public void setSwitchToOff()
	{
		/*
		 * Turn off the reactions that were previously on.
		 */
		for ( int aReac : getSpeciesParam().onStateReactions )
			switchOffreaction(allReactions[aReac]);
		/*
		 * Turn on the reactions that should now be on.
		 */
		for ( int aReac : getSpeciesParam().offStateReactions )
			switchOnReaction(allReactions[aReac]);
		/*
		 * 
		 */
		setSwitchState(false);
	}

	/**
	 * \brief Return a boolean noting whether the switch is on.
	 * 
	 * @return	Boolean stating the current switch state.
	 */
	public Boolean switchIsOn()
	{
		return switchState;
	}
	
	/**
	 * \brief Set the switch to the state in the specified boolean input.
	 * 
	 * @param newState	The state which the switch should be set to (true for
	 * on, false for off).
	 */
	public void setSwitchState(Boolean newState)
	{
		switchState = newState;
	}
	
	/**
	 * \brief Prevention of switching during the initial adaption lag period.
	 * 
	 * Returns false if this lag needs to be implemented, or true if ready to
	 * switch.
	 * 
	 * @return	Boolean noting whether the reaction switch value can be
	 * changed.
	 */
	public Boolean readyToSwitch() 
	{
		/*
		 * Don't switch if we're still within the lag.
		 * 
		 * TODO: make this a distribution around the time rather than
		 * discrete decision.
		 */
		Double stateTime = SimTimer.getCurrentTime();
		Double lagTime;
		if ( switchState )
		{
			/*
			 * If it's on already, need to consider off-state lag.
			 */
			stateTime -= timeOfRequestToSwitchOff;
			lagTime = getSpeciesParam().lagSwitchOff;
		}
		else
		{
			/*
			 * If it's on already, need to consider on-state lag.
			 */
			stateTime -= timeOfRequestToSwitchOn;
			lagTime = getSpeciesParam().lagSwitchOn;
		}
		return stateTime >= lagTime;
	}

	/**
	 * \brief Return the set of parameters that is associated with the object
	 * of this species.
	 * 
	 * @return Object of BactAdaptableParam that stores the parameters 
	 * associated with this species.
	 */
	@Override
	public BactAdaptableParam getSpeciesParam() 
	{
		return (BactAdaptableParam) _speciesParam;
	}
	
	/* ---------------------- For agent_state output ---------------------- */
	
	/**
	 * \brief Specifies the header of the columns of output information for
	 * this agent.
	 * 
	 * Used in creation of results files.
	 * 
	 * @return	String specifying the header of each column of results
	 * associated with this agent.
	 */
	@Override
	public StringBuffer sendHeader() 
	{
		return super.sendHeader().append(
				",state,turnOn,turnOff,timeRequestTurnOn,timeRequestTurnOff");
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
	public StringBuffer writeOutput()
	{
  		/*
  		 * Write the data matching the header file.
  		 */
		StringBuffer tempString = super.writeOutput().append(",");
  		/*
  		 * Switch state and timing info.
  		 */
  		tempString.append(switchState ? "1," : "0,");
  		tempString.append(turnSwitchOn ? "1," : "0,");
  		tempString.append(turnSwitchOff ? "1," : "0,");
		tempString.append(timeOfRequestToSwitchOn + "," +
										timeOfRequestToSwitchOff);
		return tempString;
	}
	
	/* ------------------------ For POV-Ray output ------------------------ */
	
	/**
	 * \brief Used in the creation of POV-Ray files, replacing the version in
	 * LocatedAgent, so that information on the switch is included.
	 * 
	 * @return	String to be included in the POV-Ray output file, consisting
	 * of species name and state of reaction switch.
	 */
	@Override
	public String getName()
	{
		return _species.speciesName + (switchState ? "_ON" : "_OFF");
	}
	
	/**
	 * \brief Used for POV-Ray output, defines the colour that this species of
	 * BactAdaptable has been assigned.
	 * 
	 * @return Color object that this species of BactAdaptable has been
	 * assigned.
	 */
	@Override
	public Color getColor()
	{
		return getSpeciesParam().getCellColor(switchState);
	}

	/**
	 * \brief Writes a colour definition to the passed-in POV-Ray output
	 * stream, ensuring switch state can be represented by a different colour.
	 *  
	 * Meant for later use in macros. This overrules the version in
	 * SpecialisedAgent.
	 * 
	 * @param fr	POV-Ray output stream where this definition should be
	 * written to.
	 * @throws IOException	Exception thrown if this stream is not open.
	 */
	@Override
	public void writePOVColorDefinition(FileWriter fr) throws IOException 
	{
		BactAdaptableParam param = getSpeciesParam();
		/*
		 * Cells with switch on.
		 */
		fr.write("#declare "+_species.speciesName+"_ON = color rgb < ");
		fr.write((param.onColor.getRed()) / 255.0 + " , ");
		fr.write((param.onColor.getGreen()) / 255.0 + " , ");
		fr.write((param.onColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
		/*
		 * Cells with switch off.
		 */
		fr.write("#declare "+_species.speciesName+"_OFF = color rgb < ");
		fr.write((param.offColor.getRed()) / 255.0 + " , ");
		fr.write((param.offColor.getGreen()) / 255.0 + " , ");
		fr.write((param.offColor.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}
}
