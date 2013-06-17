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

import java.awt.Color;
import java.io.Serializable;
import simulator.geometry.ContinuousVector;
import simulator.agent.LocatedAgent;

/**
 * \brief Used by Povray3DScene to create a capsule object in a format that can be displayed graphically in POV-Ray output
 * 
 * Used by Povray3DScene to create a capsule object in a format that can be displayed graphically in POV-Ray output
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class ParticleWithCapsule implements Serializable
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Vector that notes the location of a particular LocatedAgent agent
	 */
	private VectorProperty center;

	/**
	 * Radius of the agent core, but transformed into a size that can be represented in POV-Ray
	 */
	private double _radiusCore;
	
	/**
	 * Radius of the agent capsule, if applicable, but transformed into a size that can be represented in POV-Ray
	 */
	private double	_radiusCapsule;
	
	/**
	 * Vector representing the colour that will be used to represent the core of the agent
	 */
	private VectorProperty _colorCore;
	
	/**
	 * Vector representing the colour that will be used to represent the capsule of the agent
	 */
	private VectorProperty _colorCapsule;
	
	/**
	 * String denoting the name of the core of this agent
	 */
	private String _nameCore;
	
	/**
	 * String denoting the name of the capsule of this agent
	 */
	private String _nameCapsule;
	
	/**
	 * Boolean noting whether this agent has an EPS capsule (true) or not (false)
	 */
	private boolean _hasCapsule;
	
	/**
	 * Fraction of this agent that is active biomass
	 */
	private double _activeFrac;

	/**
	 * \brief Constructor that initialises this storage object, creating the required center and colour vector properties
	 * 
	 * Constructor that initialises this storage object, creating the required center and colour vector properties
	 */
	public ParticleWithCapsule() {
		center = new VectorProperty("");
		_colorCore = new VectorProperty("color rgb");
	}

	/**
	 * \brief Constructor that initialises this storage object with a specific LocatedAgent, creating the required center and colour vector properties
	 * 
	 * Constructor that initialises this storage object with a specific LocatedAgent, creating the required center and colour vector properties
	 */
	public ParticleWithCapsule(LocatedAgent p) 
	{
		center = new VectorProperty("");
		setCenter(p.getLocation());

		_colorCore = new VectorProperty("color rgb");
		setColorCore(p.getColor());
		setCoreRadius(p.getRadius(true));

		// bvm 27.1.2009 for using color definitions

		setNameCore(p.getName());
		setActiveFrac(p.getActiveFrac());

		_hasCapsule = p.hasEPS();
		//TODO
		// NOTE: if this is set to true, need to modify the agent.Species routine
		// that creates color definitions so that the '-capsule' colors are defined
		_hasCapsule = false;
		if (_hasCapsule) {
			_radiusCapsule = p.getRadius(true)/Povray3DScene.getScaling();
			_colorCapsule = new VectorProperty("rgbf");
			setColorCapsule(p.getColorCapsule());
			// bvm 27.1.2009 for using color definitions
			setNameCapsule(p.getSpecies().speciesName+"-capsule");
		}
	}

	/**
	 * \brief Set the colour to be used for the core of this agent from the colour specified
	 * 
	 * Set the colour to be used for the core of this agent from the colour specified
	 * 
	 * @param c	Colour to be used to represent the core of this agent
	 */
	public void setColorCore(Color c) {
		_colorCore.setValues(((float) c.getRed()) / 255,
				((float) c.getGreen()) / 255, ((float) c.getBlue()) / 255);
	}

	/**
	 * \brief Set the colour to be used for the capsule of this agent from the colour specified
	 * 
	 * Set the colour to be used for the capsule of this agent from the colour specified
	 * 
	 * @param c	Colour to be used to represent the capsule of this agent
	 */
	public void setColorCapsule(Color c) {
		float r = ColorMaps.brightenValue(((float) c.getRed()) / 255, 0.5f);
		float g = ColorMaps.brightenValue(((float) c.getGreen()) / 255, 0.5f);
		float b = ColorMaps.brightenValue(((float) c.getBlue()) / 255, 0.5f);
		_colorCapsule.setValues(r, g, b, 0.999f);
	}

	/**
	 * \brief Set the name of the core being represented to that specified
	 * 
	 * Set the name of the core being represented to that specified
	 * 
	 * @param theName	String that should be set as the name of the core
	 */
	public void setNameCore(String theName) {
		_nameCore = theName;
	}

	/**
	 * \brief Set the name of the capsule being represented to that specified
	 * 
	 * Set the name of the capsule being represented to that specified
	 * 
	 * @param theName	String that should be set as the name of the capsule
	 */
	public void setNameCapsule(String theName) {
		_nameCapsule = theName;
	}

	/**
	 * \brief Set the fraction of this agent that is active biomass
	 * 
	 * Set the fraction of this agent that is active biomass
	 * 
	 * @param activeFrac	The percentage of this agent that is active biomass
	 */
	public void setActiveFrac(double activeFrac) {
		_activeFrac = activeFrac;
	}

	/**
	 * \brief Set the centre of the output agent to a transformation of the current agent location
	 * 
	 * Set the centre of the output agent to a transformation of the current agent location
	 * 
	 * @param c	The current location of the located agent, expressed as a vector
	 */
	public void setCenter(ContinuousVector c) {
		double s = Povray3DScene.getScaling();
		center.setValues(c.x/s, c.y/s, c.z/s);
	}

	/**
	 * \brief Set the radius of the output agent to a transformation of the current agent radius
	 * 
	 * Set the radius of the output agent to a transformation of the current agent radius
	 * 
	 * @param fs	The current radius of the located agent
	 */
	public void setCoreRadius(double fs) {
		_radiusCore = fs/Povray3DScene.getScaling();
	}
	/**
	 * \brief Represents the information about this particle with capsule as a string
	 * 
	 * Represents the information about particle with capsule as a string
	 * 
	 * @return String value summarising the information stored about this particle with capsule
	 */
	public String toString() {

		// bvm 27.1.2009: modified this output to use color definitions and
		// textures rather than pigments
		String core = "sphere {\n"
			+ "\t "	+ center + "\n"
			+ "\t "	+ _radiusCore + "\n"
			+ "\t pigment { " + _nameCore + "*" + _activeFrac + " }\n"
			+ "}\n";

		if (_hasCapsule) {
			String capsule = "sphere {\n"
				+ "\t " + center + "\n"
				+ "\t " + _radiusCapsule + "\n"
				+ "\t pigment { " + _nameCapsule + "*" + _activeFrac + " }\n"
				+ "}\n";
			return core + capsule;
		}

		return core;
	}
}
