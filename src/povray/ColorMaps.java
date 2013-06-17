/**
 * \package povray
 * \brief Package of classes used to create output files that can be processed using POV-Ray software.
 * 
 * Package of classes used to create output files that can be processed using POV-Ray software. This package is part of iDynoMiCS v1.2, 
 * governed by the CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ 
 * or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  
 * "http://www.cecill.info".
 */
package povray;
import java.awt.Color;
import java.io.Serializable;


@SuppressWarnings("serial")

/**
 * \brief Class containing methods to adjust POV-Ray display colours if necessary
 * 
 * Class containing methods to adjust POV-Ray display colours if necessary. Note that from version 1.2, all the methods that were included 
 * in this class previously (such as darken colour) have been removed as these were never used in iDynoMiCS. If these are required, these
 * can be gained from downloading a copy of ColorMaps class found in iDynoMiCS version 1.1
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public abstract class ColorMaps implements Serializable
{
	/**
	 * Increase brightness of a color by factor f, keeping saturation and hue
	 * 
	 * @param c	Display colour that is to be brightened
	 * @param f	Factor by which the colour is to be brightened. Integer between 0 and 1
	 * @return Colour object that has been brightened by the given factor
	 */
	static final public Color brighten(Color c, float f) {
		//make sure f is not negative nor greater than 1
		f = (f < 0 ? 0 : f);
		f = (f > 1 ? 1.0f : f);
		float r = (float) c.getRed() / 255.0f;
		float g = (float) c.getGreen() / 255.0f;
		float b = (float) c.getBlue() / 255.0f;
		return new Color(brightenValue(r, f), brightenValue(g, f),
				brightenValue(b, f));
	}
	
	/**
	 * Creates a new section of an RGB colour by adjusting the colour value by a given factor
	 * 
	 * @param x	Section of the colur (R, G, or B) that is to be adjusted
	 * @param f	Adjustment to be applied to that section of the colour
	 * @return Float object giving the level of either R,G,or B after adjustment
	 */
	static final public float brightenValue(float x, float f) {
		return x + (1-x)*f;
	}

}
