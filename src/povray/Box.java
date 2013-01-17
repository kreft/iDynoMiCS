
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * 
 */

package povray;

import java.io.Serializable;

public class Box implements Serializable{
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private VectorProperty corner1;
	private VectorProperty corner2;
	private VectorProperty color;

	public Box() {
		corner1 = new VectorProperty("");
		corner2 = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	/**
	 * @param fs
	 */
	public void setCorner1(double x, double y, double z) {
		corner1.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setCorner2(double x, double y, double z) {
		corner2.setValues(x, y, z);
	}

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
