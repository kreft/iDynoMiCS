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

/** Parameters common to all instances of a same species */
public class EpiBacParam extends BactEPSParam {

	// number of cells contacted per hour
	public double scanSpeed;

	// parameters for growth dependency: values specify transition points
	public double lowTonusCutoff;
	public double highTonusCutoff;

	// colors for povray output
	public Color dColor, tColor, rColor;

	public EpiBacParam() {
		super();
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) {
		super.init(aSim, aSpeciesRoot,speciesDefaults);
		String colorName;

		scanSpeed = aSpeciesRoot.getParamDbl("scanSpeed");

		// if these values are not input, make them large negative to have no effect
		lowTonusCutoff = aSpeciesRoot.getParamDbl("lowTonusCutoff");
		if (Double.isNaN(lowTonusCutoff)) lowTonusCutoff = -Double.MAX_VALUE;
		highTonusCutoff = aSpeciesRoot.getParamDbl("highTonusCutoff");
		if (Double.isNaN(highTonusCutoff)) highTonusCutoff = -Double.MAX_VALUE;

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
