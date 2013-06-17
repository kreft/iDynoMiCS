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
 * \brief Class used to set the properties of the box used in POV-Ray output
 * 
 * Class used to set the properties of the box used in POV-Ray output
 *
 */
public class Box implements Serializable
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Vector property noting the coordinates of the first corner
	 */
	private VectorProperty corner1;
	
	/**
	 * Vector property noting the coordinates of the second corner
	 */
	private VectorProperty corner2;
	
	/**
	 * The colour that this box will be
	 */
	private VectorProperty color;

	/**
	 * \brief Constructor to set the coordinates and colour of this box
	 * 
	 * Constructor to set the coordinates and colour of this box
	 */
	public Box() {
		corner1 = new VectorProperty("");
		corner2 = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * \brief Set the colour that will be used to represent this box
	 * 
	 * Set the colour that will be used to represent this box
	 * 
	 * @param r	The amount of red in the colour required
	 * @param g	The amount of green in the colour required
	 * @param b	The amount of blue in the colour required
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	/**
	 * \brief Set the coordinates of the first corner of this box
	 * 
	 * Set the coordinates of the first corner this box
	 * 
	 * @param x	Double value X coordinate
	 * @param y Double value Y coordinate
	 * @param z Double value Z coordinate
	 */
	public void setCorner1(double x, double y, double z) {
		corner1.setValues(x, y, z);
	}

	/**
	 * \brief Set the coordinates of the second corner of this box
	 * 
	 * Set the coordinates of the second corner this box
	 * 
	 * @param x	Double value X coordinate
	 * @param y Double value Y coordinate
	 * @param z Double value Z coordinate
	 */
	public void setCorner2(double x, double y, double z) {
		corner2.setValues(x, y, z);
	}

	/**
	 * \brief Summarise the characteristics of this box as a string
	 * 
	 * Summarise the characteristics of this box as a string
	 * 
	 * @return String value containing all the attributes of this box
	 */
	public String toString() {
		return "box {\n"
			+ "\t "
			+ corner1
			+ "\n"
			+ "\t "
			+ corner2
			+ "\n"
			+ "\t pigment { "
			+ color
			+ " }\n"
			+ "\t\tfinish {\n"
			+ "\t\t\t phong 0.9\n"
			+ "\t\t\t phong_size 60\n"
			+ "\t\t metallic }\n"
			+ "}\n";
	}
}
