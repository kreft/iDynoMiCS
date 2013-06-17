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
 * \brief Class used to specify the light source information for POV-Ray images
 * 
 * Class used to specify the light source information for POV-Ray images
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class LightSource implements Serializable
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Vector property to hold the coordinates of the light source
	 */
	private VectorProperty location;
	
	/**
	 * Vector property to hold the details of the colour of this source, if applicable
	 */
	private VectorProperty color;

	/**
	 * \brief Constructor to create a light source and initialise the location and colour vector properties
	 * 
	 * Constructor to create a light source and initialise the location and colour vector properties
	 */
	public LightSource() {
		location = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * \brief Set the location of this light source in the image
	 * 
	 * Set the location of this light source in the image
	 * 
	 * @param x	Double value X coordinate
	 * @param y	Double value Y coordinate
	 * @param z	Double value Z coordinate
	 */
	public void setLocation(double x, double y, double z) {
		location.setValues(x, y, z);
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
	 * \brief Represents the information about the lightsource as a string
	 * 
	 * Represents the information about the lightsource as a string
	 * 
	 * @return String value summarising the information stored about this lightsource
	 */
	public String toString() {
		return "light_source {\n"
		+ "\t "
		+ location
		+ "\n\t"
		+ color
		+ "\n"
		+ "}\n";
	}
}
