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
 * \brief Creates the 'camera' settings required to produce POV-Ray output
 * 
 * Creates the 'camera' settings required to produce POV-Ray output
 *
 */
public class Camera implements Serializable
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Location of the camera in the view of the biofilm
	 */
	private VectorProperty location;
	
	/**
	 * Coordinates representing a look upwards
	 */
	private VectorProperty up;
	
	/**
	 * Coordinates representing a look right
	 */
	private VectorProperty right;
	
	/**
	 * Coordinates that should be focused upon
	 */
	private VectorProperty look_at;
	
	/**
	 * Angle to use to view the image
	 */
	private double angle;

	/**
	 * \brief Constructor used to initialise each of the location, up, right, and look at vector properties
	 * 
	 * Constructor used to initialise each of the location, up, right, and look at vector properties
	 */
	public Camera() 
	{
		location = new VectorProperty("location");
		up = new VectorProperty("up");
		right = new VectorProperty("right");
		look_at = new VectorProperty("look_at");
	}

	/**
	 * \brief Set the location of the camera to the specified X,Y, and Z coordinates
	 * 
	 * Set the location of the camera to the specified X,Y, and Z coordinates
	 * 
	 * @param x	Double value X coordinate
	 * @param y	Double value Y coordinate
	 * @param z	Double value Z coordinate
	 */
	public void setLocation(double x, double y, double z) {
		location.setValues(x, y, z);
	}

	/**
	 * \brief Set the upward view of the camera to the specified X,Y, and Z coordinates
	 * 
	 * Set the upward view of the camera to the specified X,Y, and Z coordinates
	 * 
	 * @param x	Double value X coordinate
	 * @param y	Double value Y coordinate
	 * @param z	Double value Z coordinate
	 */
	public void setUp(double x, double y, double z) {
		up.setValues(x, y, z);
	}

	/**
	 * \brief Set the right view of the camera to the specified X,Y, and Z coordinates
	 * 
	 * Set the right view of the camera to the specified X,Y, and Z coordinates
	 * 
	 * @param x	Double value X coordinate
	 * @param y	Double value Y coordinate
	 * @param z	Double value Z coordinate
	 */
	public void setRight(double x, double y, double z) {
		right.setValues(x, y, z);
	}

	/**
	 * \brief Set the focus of the camera to the specified X,Y, and Z coordinates
	 * 
	 * Set the focus of the camera to the specified X,Y, and Z coordinates
	 * 
	 * @param x	Double value X coordinate
	 * @param y	Double value Y coordinate
	 * @param z	Double value Z coordinate
	 */
	public void setLook_at(double x, double y, double z) {
		look_at.setValues(x, y, z);
	}

	/**
	 * \brief Set the angle of view of the camera
	 * 
	 * Set the angle of view of the camera
	 * 
	 * @param f	Double value angle of the camera view
	 */
	public void setAngle(double f) {
		angle = f;
	}

	/**
	 * \brief Represents the information about the camera as a string
	 * 
	 * Represents the information about the camera as a string
	 * 
	 * @return String value summarising the information stored about this camera
	 */
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
