/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * Bulk is an object used to define the environment connec
 * Implements 3D vector of continuous spatial coordinates
 * Can be used to store Continuous coordinates or Movement vectors
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 */

package simulator.geometry;

import java.io.Serializable;

import utils.ExtraMath;
import org.jdom.*;

public class ContinuousVector implements Cloneable, Serializable {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	public double             x, y, z;

	/* _____________________ CONSTRUCTOR _____________________________ */
	/**
	 * Default Constructor (null coordinates)
	 */
	public ContinuousVector() {
		x = 0;
		y = 0;
		z = 0;
	}

	public ContinuousVector(ContinuousVector aCC) {
		x = aCC.x;
		y = aCC.y;
		z = aCC.z;
	}

	public ContinuousVector(Element xmlRoot) {
		x = Double.parseDouble(xmlRoot.getAttributeValue("x"));
		y = Double.parseDouble(xmlRoot.getAttributeValue("y"));
		z = Double.parseDouble(xmlRoot.getAttributeValue("z"));
	}

	/**
	 * Translate a discrete coordinates expressed on a discrete spatial grid
	 * with the resolution res to a continuous one
	 * @param cc
	 * @param res
	 */
	public ContinuousVector(DiscreteVector dC, double res) {
		x = (.5+dC.i)*res;
		y = (.5+dC.j)*res;
		z = (.5+dC.k)*res;
	}

	/**
	 * Parametrized constructor.
	 * @param x
	 * @param y
	 * @param z
	 */
	public ContinuousVector(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public void set(ContinuousVector cc) {
		this.x = cc.x;
		this.y = cc.y;
		this.z = cc.z;
	}

	public void set(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/* ___________________________ TOOLS ___________________________________ */
	/**
	 * set x = y = z = 0
	 */
	public void reset() {
		x = 0;
		y = 0;
		z = 0;
	}

	public boolean isValid() {
		return !(Double.isInfinite(x)|Double.isNaN(x)|Double.isInfinite(y)|Double.isNaN(y)
		        |Double.isInfinite(z)|Double.isNaN(z));
	}

	/**
	 * Changes the sign of the vector. Used for movement vectors.
	 */
	public void turnAround() {
		x = -x;
		y = -y;
		z = -z;
	}

	/**
	 * @return identity to coordinate (x,y,z). true = is identical
	 */
	public boolean equals(double x, double y, double z) {
		if (this.x==x&this.y==y&this.z==z) return true;
		else return false;
	}

	public boolean isZero() {
		return (this.x==0&this.y==0&&this.z==0);
	}

	/**
	 * Print coordinates to string
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return ExtraMath.toString(x, false)+",\t"+ExtraMath.toString(y, false)+",\t"
		        +ExtraMath.toString(z, false);
	}

	/**
	 * add vector v
	 */
	public void add(ContinuousVector v) {
		this.x += v.x;
		this.y += v.y;
		this.z += v.z;
	}

	public void add(double x, double y, double z) {
		this.x += x;
		this.y += y;
		this.z += z;
	}

	public void sendSum(ContinuousVector a, ContinuousVector b) {
		this.x = a.x+b.x;
		this.y = a.y+b.y;
		this.z = a.z+b.z;
	}

	/**
	 * subtract vector v
	 */
	public void subtract(ContinuousVector v) {
		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;
	}

	/**
	 * Set to the current coordinates the difference a-b
	 */
	public void sendDiff(ContinuousVector a, ContinuousVector b) {
		this.x = a.x-b.x;
		this.y = a.y-b.y;
		this.z = a.z-b.z;
	}

	/**
	 * @return scalar product (dot product) with vector cc
	 */
	public double prodScalar(ContinuousVector cc) {
		return this.x*cc.x+this.y*cc.y+this.z*cc.z;
	}

	/**
	 * multiply (stretch) by multiplier
	 */
	public void times(double multiplier) {
		this.x *= multiplier;
		this.y *= multiplier;
		this.z *= multiplier;
	}

	/**
	 * Set to random coordinate (x,y,z)
	 */
	public void alea() {
		this.x = ExtraMath.getUniRand();
		this.y = ExtraMath.getUniRand();
		this.z = ExtraMath.getUniRand();
	}

	public void alea(boolean is3D) {
		this.x = ExtraMath.getUniRand(-1,1);
		this.y = ExtraMath.getUniRand(-1,1);
		this.z = (is3D ? ExtraMath.getUniRand(-1,1) : 0);
	}

	/**
	 * Normalize this Vector to unit length.
	 * 
	 * @param newLength
	 */
	public void normalizeVector() {
		double v = this.norm();
		if (v!=0) this.times(1/this.norm());

	}

	/**
	 * Normalize this Vector to a given length.
	 * 
	 * @param newLength
	 */
	public void normalizeVector(double newLength) {
		this.times(newLength/this.norm());
	}

	/**
	 * @return absolute distance to cc
	 */
	public double distance(ContinuousVector cc) {
		return (double) Math.sqrt(Math.abs((this.x-cc.x))*Math.abs((this.x-cc.x))+Math.abs((this.y-cc.y))*Math.abs((this.y-cc.y))
		        + Math.abs((this.z-cc.z))*Math.abs((this.z-cc.z)));
	}

	/**
	 * @return norm (absolute length)
	 */
	public double norm() {
		return Math.sqrt(x*x+y*y+z*z);
	}

	/**
	 * @return cosine of the angle to vector cc
	 */
	public double cosAngle(ContinuousVector cc) {
		return (x*cc.x+y*cc.y+z*cc.z)/Math.sqrt((x*x+y*y+z*z)*(cc.x*cc.x+cc.y*cc.y+cc.z*cc.z));
	}

	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

}