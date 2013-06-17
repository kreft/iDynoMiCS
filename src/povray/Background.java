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

import java.io.Serializable;

/**
 * \brief Specifies the background that will be used in the POV-Ray images
 * 
 * Specifies the background that will be used in the POV-Ray images
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class Background implements Serializable
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Colour that this background will be assigned, specified as a vector of red, green, and blue colours
	 */
	private VectorProperty color;

	/**
	 * \brief Constructor to set the colour as a new RGB vector property
	 * 
	 * Constructor to set the colour as a new RGB vector property
	 */
	public Background() {
		color = new VectorProperty("color rgb");
	}

	/**
	 * \brief Sets the value of this vector property to reflect a colour, specified in red, green, and blue values
	 * 
	 * Sets the value of this vector property to reflect a colour, specified in red, green, and blue values
	 * 
	 * @param r	The amount of red in the colour required
	 * @param g	The amount of green in the colour required
	 * @param b	The amount of blue in the colour required
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	/**
	 * \brief Return the vector property as a string denoting what colour the vector represents
	 * 
	 * Return the vector property as a string denoting what colour the vector represents
	 * 
	 * @return String value denoting what colour the vector represents
	 */
	public String toString() {
		return "background {\n"
		+ "\t" + color + "\n"
		+ "}\n";
	}
}
