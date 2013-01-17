/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * The scene background
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

public class Background implements Serializable{
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private VectorProperty color;

	public Background() {
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	public String toString() {
		return "background {\n"
		+ "\t" + color + "\n"
		+ "}\n";
	}
}
