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

import simulator.Simulator;



import java.awt.Color;

import utils.XMLParser;

/**
 * Parameters common to all instances of a same species
 */
public class EpiBacParam extends BactEPSParam
{
	/**
	 * Number of cells contacted (per hour).
	 */
	public Double scanSpeed;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public Double lowTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public Double highTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Colours for POV-Ray output.
	 */
	public Color dColor, tColor, rColor;
	
	public EpiBacParam()
	{
		super();
	}
	
	/**
	 * 
	 */
	@Override
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot,speciesDefaults);
		
		/*
		 * Get the scanSpeed and two tonus parameters.
		 */
		Double temp;
		
		temp = getSpeciesParameterDouble("scanSpeed", aSpeciesRoot,
															speciesDefaults);
		scanSpeed = Double.isFinite(temp) ? temp : scanSpeed;
		
		temp = getSpeciesParameterDouble("lowTonusCutoff", aSpeciesRoot,
															speciesDefaults);
		lowTonusCutoff = Double.isFinite(temp) ? temp : lowTonusCutoff;
		
		temp = getSpeciesParameterDouble("highTonusCutoff", aSpeciesRoot,
															speciesDefaults);
		highTonusCutoff = Double.isFinite(temp) ? temp : highTonusCutoff;
		
		/*
		 * Get the colours for POV-Ray output.
		 */
		String colorName;
		
		colorName = aSpeciesRoot.getParam("donorColor");
		if (colorName == null)
			colorName = "white";
		dColor = utils.UnitConverter.getColor(colorName);

		colorName = aSpeciesRoot.getParam("transconjugantColor");
		if (colorName == null)
			colorName = "white";
		tColor = utils.UnitConverter.getColor(colorName);

		colorName = aSpeciesRoot.getParam("recipientColor");
		if (colorName == null)
			colorName = "white";
		rColor = utils.UnitConverter.getColor(colorName);
	}
}
