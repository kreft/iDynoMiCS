
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Vector property for PovRay descriptive language
 * 
 */

/**
 * @since Feb 2007
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 */

package povray;

import java.io.Serializable;


public class VectorProperty  implements Serializable {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private String _name;
	private double[] _values;

	public VectorProperty(String name) {
		_name = name;
	}

	public void setValues(double x, double y, double z) {
		_values = new double[3];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
	}

	public void setValues(double x, double y, double z, double w) {
		_values = new double[4];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
		_values[3] = w;
	}

	public String toString() {
		int n = _values.length;
		String out = _name + " <";
		for (int i = 0; i < n; i++) {
			out += " " + _values[i] + (i == n-1 ? " >" : ", ");
		}
		return out;
	}
}
