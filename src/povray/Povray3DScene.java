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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import simulator.Simulator;
import simulator.agent.Species;
import simulator.geometry.Domain;

/**
 * \brief Create a POV-Ray 3D scene object for use in visualising 3D Biofilms (with class Biofilm3D)
 * 
 * Create a POV-Ray 3D scene object for use in visualising 3D Biofilms (with class Biofilm3D)
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class Povray3DScene implements Serializable 
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long   serialVersionUID = 1L;

	/**
	 * Current simulator object through which the result has been formed
	 */
	public Simulator            mySim;
	
	final private static String INCHEADER        = "sceneheader.inc";
	final private static String INCFOOTER        = "scenefooter.inc";

	/**
	 * Camera object specifying the 'camera' settings required to produce POV-Ray output
	 */
	private Camera              _camera;
	
	/**
	 * Specifies the background that will be used in the POV-Ray images
	 */
	private Background          _background;
	
	/**
	 * Used to specify the light source information for POV-Ray images
	 */
	private LightSource[]       _lightSource;
	
	/**
	 * Biofilm object that contains the boundary conditions and cells that will be drawn on the output
	 */
	private Biofilm3D           _biofilm;
	
	/**
	 * The simulation computation domain that will be represented
	 */
	private Domain              _domain;

	/**
	 * Scaling that needs to be applied to represent the computation domain on the output
	 */
	private static double       _scaling;
	
	/**
	 * Scaled X value used to scale other X coordinates that are to be included on the output
	 */
	private double              _x;
	
	/**
	 * Scaled Y value used to scale other X coordinates that are to be included on the output
	 */
	private double				_y;
	
	/**
	 * Scaled Z value used to scale other X coordinates that are to be included on the output
	 */
	private double				_z;
	
	/**
	 * \brief Constructor that initialises the POV-Ray scene, setting the domain to be represented and the rendering scale that will be used
	 * 
	 * Constructor that initialises the POV-Ray scene, setting the domain to be represented and the rendering scale that will be used
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param domainName	The name of the computation domain to be rendered
	 */
	public Povray3DScene(Simulator aSim,String domainName) {
		mySim = aSim;
		_domain = aSim.world.getDomain(domainName);
		_scaling = _domain.getLongestSize();
		_x = _domain.length_X/_scaling;
		_y = _domain.length_Y/_scaling;
		_z = _domain.length_Z/_scaling;

		initializeScene();
	}

	/**
	 * \brief Initialise scene parameters (camera, background and light)
	 * 
	 * Initialise scene parameters (camera, background and light)
	 */
	private void initializeScene() {
		// Set the camera
		_camera = new Camera();
		_lightSource = new LightSource[2];

		// location and orientation (where is and where it looks)
		// NOTE THAT POVRAY USES A LEFT-HAND COORDINATE SYSTEM
		// (x is left/right, y is up/down, and z is into screen;
		//  this means compared to the idyno representation, povray swaps x and y)

		if (_domain.is3D){
			_camera.setLocation(0, _y*1.5, -_x*1.5);

			_lightSource[0] = new LightSource();
			_lightSource[0].setLocation(_y, _x, -_z);
			_lightSource[0].setColor(1f, 1f, 1f);

			_lightSource[1] = new LightSource();
			_lightSource[1].setLocation(-_y, _x, -_z);	
			_lightSource[1].setColor(1f, 1f, 1f);		
		}else{
			_camera.setLocation(0, 0, -_x*1.5);

			_lightSource[0] = new LightSource();
			_lightSource[0].setLocation(_z, 0, -_x*1.5);
			_lightSource[0].setColor(1f, 1f, 1f);

			_lightSource[1] = new LightSource();
			_lightSource[1].setLocation(-_z, 0, -_x*1.5);
			_lightSource[1].setColor(1f, 1f, 1f);
		}

		_camera.setLook_at(0, 0, 0);

		// these set the aspect ratio of the view
		_camera.setUp(0, 1, 0);
		_camera.setRight(-1.33f, 0, 0);
		_camera.setAngle(60);

		// Set the background
		_background = new Background();
		_background.setColor(1f, 1f, 1f);

	}


	/**
	 * \brief Write the complete model state to a file (including camera, lights and boxes)
	 * 
	 * Write the complete model state to a file (including camera, lights and boxes)
	 * 
	 * @param fileName	File name of the POV-Ray file being written
	 * @throws IOException	Exception thrown if the file output stream cannot be opened
	 */
	public void writeModelStateFull(String fileName) throws IOException 
	{

		File sysData = new File(fileName);
		FileWriter fr = new FileWriter(sysData);

		fr.write(_camera.toString());
		fr.write(_background.toString());
		fr.write(_lightSource[0].toString());
		fr.write(_lightSource[1].toString());

		_biofilm = new Biofilm3D(this);
		_biofilm.modelStateToFile(fr);

		fr.close();
	}

	/**
	 * \brief Write the include files. For use with writeModelState
	 * 
	 * Write the include files. For use with writeModelState
	 * 
	 * @param dir Results directory to write the include file to
	 * @throws IOException	Exception thrown if the file output stream cannot be opened	
	 */
	public File[] writePovrayIncFiles(String dir) throws IOException 
	{
		_biofilm = new Biofilm3D(this);

		// header include file
		java.io.File header = new java.io.File(dir+INCHEADER);
		FileWriter fr = new FileWriter(header);
		fr.write(_camera.toString());
		fr.write(_background.toString());
		fr.write(_lightSource[0].toString());
		fr.write(_lightSource[1].toString());
		_biofilm.biofilmHeaderToFile(fr);

		// bvm 27.1.2009: define colors based on species name; for use in macros
		// (uses the progenitor, which is able to write colors for different individual states)
		for (Species aSpecies : mySim.speciesList)
			aSpecies.getProgenitor().writePOVColorDefinition(fr);

		fr.close();

		// footer include file
		java.io.File footer = new java.io.File(dir+INCFOOTER);
		fr = new FileWriter(footer);
		_biofilm.biofilmFooterToFile(fr);
		fr.close();

		// prepare file array to return
		File[] incs = { header, footer };
		return incs;
	}

	/**
	 * \brief Write the present state using include files for camera, background lights and solid surface. 
	 * 
	 * Write the present state using include files for camera, background lights and solid surface. Using this method instead of 
	 * writeModelStateFull, which includes the full scene information in each .pov file created at an iteration, allows changing scene 
	 * properties after the simulation
	 * 
	 * @param fileName	Filename of the POV-Ray file to be written
	 * @return the POV-Ray file containing the model state
	 * @throws IOException	Exception thrown should this file not be able to be accessed
	 */
	public File writeModelState(String fileName) throws IOException 
	{
		_biofilm = new Biofilm3D(this);

		File sysData = new File(fileName);
		FileWriter fr = new FileWriter(sysData);


		fr.write("#include \""+"colors.inc"+"\"\n");

		fr.write("#include \""+INCHEADER+"\"\n");

		try {
			_biofilm.particlesToFile(fr);
		} catch(Exception e) {
			_biofilm.particlesToFile(fr);
		}

		fr.write("#include \""+INCFOOTER+"\"\n");
		fr.close();
		return sysData;
	}

	/**
	 * \brief Make top perspective
	 * 
	 * Make top perspective
	 */
	public void setTopPerspective() {
		_camera.setLocation(0, _y*1.56f, 0);
		_lightSource[1].setLocation(_z, 0, 0);
	}

	/**
	 * \brief Make side perspective
	 * 
	 * Make side perspective
	 */
	public void setSidePerspective() {
		_camera.setLocation(0, _y, _x*1.17f);
		_lightSource[1].setLocation(_z, 0, 0);
		_camera.setLook_at(0, -_y*0.5, 0);
	}

	/**
	 * \brief Make angle perspective
	 * 
	 * Make angle perspective
	 */
	public void setAnglePerspective() {
		_camera.setLocation(_z, _y, _x*2.5f);
		_lightSource[1].setLocation(_z, 0, 0);
		_camera.setLook_at(0, 0, 0);
	}

	/**
	 * \brief Add a cell to the biofilm, to be used with scenes constructed using the Povray3DScene(float x, float y, float z, int n) constructor
	 * 
	 * Add a cell to the biofilm, to be used with scenes constructed using the Povray3DScene(float x, float y, float z, int n) 
	 * constructor
	 * 
	 * @param x	X coordinate of the cell to add
	 * @param y	Y coordinate of the cell to add
	 * @param z	Z coordinate of the cell to add
	 * @param rad	Radius of the cell to add
	 * @param r	Amount of red colour to use in creating the colour of this cell
	 * @param g	Amount of green colour to use in creating the colour of this cell
	 * @param b	Amount of blue colour to use in creating the colour of this cell
	 */
	public void addCell(float x, float y, float z, float rad, int r, int g, int b) {
		_biofilm.addCell(x, y, z, rad, r, g, b);
	}

	/**
	 * \brief Get the scaled X value used to scale other X coordinates
	 * 
	 * Get the scaled X value used to scale other X coordinates
	 * 
	 * @return Scaled X value
	 */
	protected double getX() {
		return _x;
	}

	/**
	 * \brief Get the scaled Y value used to scale other Y coordinates
	 * 
	 * Get the scaled Y value used to scale other Y coordinates
	 * 
	 * @return Scaled Y value
	 */
	protected double getY() {
		return _y;
	}

	/**
	 * \brief Get the scaled Z value used to scale other Z coordinates
	 * 
	 * Get the scaled Z value used to scale other Z coordinates
	 * 
	 * @return Scaled Z value
	 */
	protected double getZ() {
		return _z;
	}

	/**
	 * \brief Return the computational domain being represented in this POV-Ray output
	 * 
	 * Return the computational domain being represented in this POV-Ray output
	 * 
	 * @return	Computation domain being represented in this POV-Ray output
	 */
	Domain getDomain() 
	{
		return _domain;
	}

	/**
	 * \brief Return the scaling factor applied to points in this output
	 * 
	 * Return the scaling factor applied to points in this output
	 * 
	 * @return Calculated scaling value
	 */
	public static double getScaling() {
		return _scaling;
	}

	/**
	 * @param os
	 * @throws IOException
	 */
	//public static void serializeStaticState(ObjectOutputStream os) throws IOException {
		//os.writeDouble(_scaling);
	//}

	/**
	 * @param os
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	//public static void deserializeStaticState(ObjectInputStream os) throws IOException,
	//ClassNotFoundException {
		//_scaling = os.readFloat();
	//}

}
