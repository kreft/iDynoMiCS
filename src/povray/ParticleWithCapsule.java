
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

/**
 * @since Feb 2007
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */


package povray;

import java.awt.Color;
import java.io.Serializable;

import simulator.geometry.ContinuousVector;
import simulator.agent.LocatedAgent;


public class ParticleWithCapsule implements Serializable{
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private VectorProperty center;

	private double _radiusCore, _radiusCapsule;
	private VectorProperty _colorCore,_colorCapsule;
	private String _nameCore, _nameCapsule;
	private boolean _hasCapsule;
	private double _activeFrac;

	/* _________________ CONSTRUCTOR _________________________ */
	public ParticleWithCapsule() {
		center = new VectorProperty("");
		_colorCore = new VectorProperty("color rgb");
	}

	public ParticleWithCapsule(LocatedAgent p) {
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
	 * @param color
	 */
	public void setColorCore(Color c) {
		_colorCore.setValues(((float) c.getRed()) / 255,
				((float) c.getGreen()) / 255, ((float) c.getBlue()) / 255);
	}

	/**
	 * For now sets capsule to gray
	 * 
	 * @param fs
	 */
	public void setColorCapsule(Color c) {
		float r = ColorMaps.brightenValue(((float) c.getRed()) / 255, 0.5f);
		float g = ColorMaps.brightenValue(((float) c.getGreen()) / 255, 0.5f);
		float b = ColorMaps.brightenValue(((float) c.getBlue()) / 255, 0.5f);
		_colorCapsule.setValues(r, g, b, 0.999f);
	}

	/**
	 * @param theName
	 */
	public void setNameCore(String theName) {
		_nameCore = theName;
	}

	/**
	 * @param theName
	 */
	public void setNameCapsule(String theName) {
		_nameCapsule = theName;
	}

	/**
	 * @param activeFrac
	 */
	public void setActiveFrac(double activeFrac) {
		_activeFrac = activeFrac;
	}

	/**
	 * @param fs
	 */
	public void setCenter(ContinuousVector c) {
		double s = Povray3DScene.getScaling();
		center.setValues(c.x/s, c.y/s, c.z/s);
	}

	/**
	 * @param fs
	 */
	public void setCoreRadius(double fs) {
		_radiusCore = fs/Povray3DScene.getScaling();
	}

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
