
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

/**
 * @since Feb 2007
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */


package povray;

import java.io.Serializable;




public class LightSource implements Serializable{
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private VectorProperty location;
	private VectorProperty color;

	/**
	 * 
	 */
	public LightSource() {
		location = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setLocation(double x, double y, double z) {
		location.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

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
