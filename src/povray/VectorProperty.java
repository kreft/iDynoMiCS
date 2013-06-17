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
 * \brief Class used by the majority of POV-Ray scripts to represent coordinates and colours as vector properties for use in POV-Ray
 * 
 * Class used by the majority of POV-Ray scripts to represent coordinates and colours as vector properties for use in POV-Ray
 *
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 */
public class VectorProperty  implements Serializable 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Name assigned to this vector property
	 */
	private String _name;
	
	/**
	 * Double array value that this object will take
	 */
	private double[] _values;

	/**
	 * \brief Constructor to create the object and give it a specific name
	 * 
	 * Constructor to create the object and give it a specific name
	 * 
	 * @param name	Name to be assigned to this vector property object
	 */
	public VectorProperty(String name) {
		_name = name;
	}

	/**
	 * \brief Set the values of this vector property, storing the values provided into the array
	 * 
	 * Set the values of this vector property, storing the values provided into the array
	 * 
	 * @param x	Either an X coordinate or an amount of red in a colour
	 * @param y Either an Y coordinate or an amount of green in a colour
	 * @param z Either an Z coordinate or an amount of blue in a colour
	 */
	public void setValues(double x, double y, double z) {
		_values = new double[3];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
	}

	/**
	 * \brief Set the values of this vector property for a object of 4 values, storing the values provided into the array
	 * 
	 * Set the values of this vector property for a object of 4 values, storing the values provided into the array
	 * 
	 * @param x	Value of this vector object
	 * @param y Value of this vector object
	 * @param z Value of this vector object
	 * @param w Value of this vector object
	 */
	public void setValues(double x, double y, double z, double w) {
		_values = new double[4];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
		_values[3] = w;
	}

	/**
	 * \brief Summarise this vector object as a string
	 * 
	 * Summarise this vector object as a string
	 * 
	 * @return String value summarising this vector object
	 */
	public String toString() {
		int n = _values.length;
		String out = _name + " <";
		for (int i = 0; i < n; i++) {
			out += " " + _values[i] + (i == n-1 ? " >" : ", ");
		}
		return out;
	}
}
