/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent.zoo;

import java.util.ArrayList;

import org.jdom.Element;
import java.awt.Color;

import simulator.Simulator;
import simulator.SoluteGrid;
import utils.XMLParser;

/**
 * \brief Creates a parameter object for BactAdaptable, a Bacterium agent object that can change its reaction based on local conditions
 * 
 * Creates a parameter object for BactAdaptable, a Bacterium agent object that can change its reaction based on local conditions
 * 
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 *
 */
public class BactAdaptableParam extends BactEPSParam 
{

	/**
	 * Array used to read the solute concentrations
	 */
	protected SoluteGrid[] _soluteList;

	
	/**
	 * Subsets of all known reactions which take place when the switch is set to on. Contains indices to be used in referencing the allReactions array
	 */
	protected ArrayList<Integer> onStateReactions;
	
	/**
	 * Subsets of all known reactions which take place when the switch is set to off. Contains indices to be used in referencing the allReactions array
	 */
	protected ArrayList<Integer> offStateReactions;

	/**
	 * The factor controlling what changes the value of this switch - the local solute concentration or agents biomass amount
	 */
	protected String switchType;
	
	/**
	 * The name of the component that is controlling the value of this switch (e.g. o2d) 
	 */
	protected String switchControl;
	
	/**
	 * Integer index of the particle or reaction in the simulation dictionary
	 */
	protected int switchControlIndex;
	
	/**
	 * Used in the onCondtion reaction mark-up - tests switch changing thresholds. Can be lessThan or greaterThan
	 */
	protected String switchCondition; 
	
	/**
	 * The value (concentration or mass) that will prompt the changing of the reaction switch
	 */
	protected double switchValue;
	
	/**
	 * Lag time before switch actually changes after a request to change it, to capture microbial metabolic delays
	 */
	protected double lagSwitchOn  = 1.0;
	
	/**
	 * Lag time before switch actually changes after a request to change it, to capture microbial metabolic delays
	 */
	protected double lagSwitchOff = 1.0;

	/**
	 * Colour to be used by POV-Ray to show objects of this species when the reaction switch is on
	 */
	public Color onColor;
	
	/**
	 * Colour to be used by POV-Ray to show objects of this species when the reaction switch is off
	 */
	public Color offColor;

	/**
	 * \brief Creates a parameter storage object for the BactAdaptable species type
	 * 
	 * Creates a parameter storage object for the BactAdaptable species type
	 */
	public BactAdaptableParam() 
	{
		super();
	}

	/**
	 * \brief Initialises BactAdaptable Species parameters, calling the relevant superclasses to initialise common parameters and then setting up the reaction switch
	 * 
	 * Initialises Bacterium EPS Species parameters, calling the relevant superclasses to initialise common parameters and then setting up the reaction switch
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		String colorName;

		super.init(aSim,aSpeciesRoot,speciesDefaults);

		// Now take care of reaction switch business
		int reacIndex;
		_soluteList = aSim.soluteList;

		onStateReactions = new ArrayList<Integer>();
		offStateReactions = new ArrayList<Integer>();

		XMLParser switchParser = new XMLParser(aSpeciesRoot.getChildElement("reactionSwitch"));
		XMLParser childParser;

		//////////////////////////////////////////////////////////////////////////
		// create list of reactions during the On state
		childParser = new XMLParser(switchParser.getChildElement("whenOn"));
		for (Element aReactionMarkUp : childParser.buildSetMarkUp("reaction")) {
			// Add the reaction to the list of reactions for the on-state
			// but add only ACTIVE reactions!!
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			if (aReactionMarkUp.getAttributeValue("status").equals("active")) {
				onStateReactions.add(reacIndex);

			}
		}

		colorName = childParser.getParam("color");
		if (colorName == null)
			colorName = "white";
		onColor = utils.UnitConverter.getColor(colorName);

		lagSwitchOn  = childParser.getParamTime("switchLag");


		//////////////////////////////////////////////////////////////////////////
		// create list of reactions during the Off state
		childParser = new XMLParser(switchParser.getChildElement("whenOff"));
		for (Element aReactionMarkUp : childParser.buildSetMarkUp("reaction")) {
			// Add the reaction to the list of reactions for the off-state
			// but add only ACTIVE reactions!!
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			if (aReactionMarkUp.getAttributeValue("status").equals("active")) {
				offStateReactions.add(reacIndex);

			}
		}

		colorName = childParser.getParam("color");
		if (colorName == null)
			colorName = "black";
		offColor = utils.UnitConverter.getColor(colorName);

		lagSwitchOff = childParser.getParamTime("switchLag");



		// now set up when the on and off states are set/met
		childParser = new XMLParser(switchParser.getChildElement("onCondition"));

		// get type (solute or biomass) and controlling component
		switchType = childParser.getAttribute("type");
		switchControl = childParser.getAttribute("name");

		// get whether it's lessThan or greaterThan
		switchCondition = childParser.getParam("switch");

		// get the concentration or mass for the switch
		if (switchType.equals("solute")) {
			switchValue = childParser.getParamConc("concentration");
			switchControlIndex = aSim.soluteDic.indexOf(switchControl);
			if (switchControlIndex == -1)
				System.out.println("WARNING: solute "+switchControl+
				" used in the <reactionSwitch> markup does not exist.");
		}
		if (switchType.equals("biomass")) {
			switchValue = childParser.getParamMass("mass");
			switchControlIndex = aSim.particleDic.indexOf(switchControl);
			if (switchControlIndex == -1)
				System.out.println("WARNING: biomass type "+switchControl+
				" used in the <reactionSwitch> markup does not exist.");

		}


	}

}
