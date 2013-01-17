
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * 
 */

package povray;

import java.io.Serializable;


public class Camera implements Serializable{
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private VectorProperty location;
	private VectorProperty up;
	private VectorProperty right;
	private VectorProperty look_at;
	private double angle;

	public Camera() {
		location = new VectorProperty("location");
		up = new VectorProperty("up");
		right = new VectorProperty("right");
		look_at = new VectorProperty("look_at");
	}

	/**
	 * @param fs
	 */
	public void setLocation(double x, double y, double z) {
		location.setValues(x, y, z);
	}

	public void setUp(double x, double y, double z) {
		up.setValues(x, y, z);
	}

	public void setRight(double x, double y, double z) {
		right.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setLook_at(double x, double y, double z) {
		look_at.setValues(x, y, z);
	}

	/**
	 * @param f
	 */
	public void setAngle(double f) {
		angle = f;
	}

	public String toString() {
		return "camera {\n"+ "\t"+ location
			+ "\n"+ "\t "
			+ up+ "\n"+ "\t "
			+ right+ "\n"+ "\t "
			+ look_at+ "\n"+ "\tangle "
			+ angle+ "\n"
			+ "}\n";
	}

}
