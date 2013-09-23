/**
 * Project iDynoMicS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Simulator is the top-level class of the simulation core.
 * It is used to create and run a simulation; this class is called by the class
 * MicroCoSim.
 * 
 * @since June 2006
 * @version 0.1
 * @author Andreas Doetsch, University of Bonn
 * @author Laurent Lardon, DTU Denmark
 */

package simulator.agent.zoo;

import simulator.Simulator;
import java.awt.Color;
import java.util.*;

import org.jdom.Element;

import utils.*;

/** Parameters common to all instances of a same species */
	public class MultiEpiBacParam extends BactEPSParam {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID     = 1L;

	public double             donorProbability    = 1;
	public double             recipientProbability = 1;

	//Plasmid replication speed
	public double             nBind                = 2;
	public double             KSat                 = 0.132;
	
	// number of cells touched by hour
	//sonia 11.10.2010 scan speed is now a plasmid parameter
	//public double             scanSpeed;
	public Color              dColor, tColor, rColor;
	
	// sonia: useful..
	public String epiBacName;
	

	
	public MultiEpiBacParam() {
		super();
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) {
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		double value;

		value = aSpeciesRoot.getParamDbl("donorProbability");
		donorProbability = (Double.isNaN(value) ? donorProbability : value);

		value = aSpeciesRoot.getParamDbl("recipientProbability");
		recipientProbability = (Double.isNaN(value) ? recipientProbability : value);

		//scanSpeed = aSpeciesRoot.getParamDbl("scanSpeed");

		String colorName;
		colorName = aSpeciesRoot.getParam("donorColor");
		if (colorName==null) colorName = "white";
		dColor = utils.UnitConverter.getColor(colorName);

		colorName = aSpeciesRoot.getParam("transconjugantColor");
		if (colorName==null) colorName = "white";
		tColor = utils.UnitConverter.getColor(colorName);

		colorName = aSpeciesRoot.getParam("recipientColor");
		if (colorName==null) colorName = "white";
		rColor = utils.UnitConverter.getColor(colorName);

		value = aSpeciesRoot.getParamDbl("nBind");
		nBind = (Double.isNaN(value) ? nBind : value);

		value = aSpeciesRoot.getParamDbl("KSat");
		KSat = (Double.isNaN(value) ? KSat : value);
		
		
		//retrieving bacteria names
		epiBacName = aSpeciesRoot.getAttribute("name");

		
	}

}
