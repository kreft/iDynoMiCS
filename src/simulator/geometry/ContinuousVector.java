/**
 * \package simulator.geometry
 * \brief Package of boundary utilities that aid the creation of the environment being simulated
 * 
 * Package of boundary utilities that aid the creation of the environment being simulated. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry;

import java.io.Serializable;
import org.jdom.*;
import utils.ExtraMath;

/**
 * \brief Implements 3D vector of continuous spatial coordinates. Can be used to store Continuous coordinates or Movement vectors
 * 
 * Implements 3D vector of continuous spatial coordinates. Can be used to store Continuous coordinates or Movement vectors
 * 
 * @author Andreas D�tsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Jo�o Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public class ContinuousVector implements Cloneable, Serializable, Comparable<ContinuousVector>
{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * X coordinate of the point contained in this vector
	 */
	public double             x;
	
	/**
	 * Y coordinate of the point contained in this vector
	 */
	public double 			  y; 
	
	/**
	 * Z coordinate of the point contained in this vector
	 */
	public double 			  z;

	
	/**
	 * \brief Default Constructor - constructs a ContinuousVector with points at 0,0,0
	 * 
	 * Default Constructor - constructs a ContinuousVector with points at 0,0,0
	 */
	public ContinuousVector() {
		x = 0;
		y = 0;
		z = 0;
	}

	/**
	 * \brief Constructs a continuous vector with points specified by a provided continuous vector
	 * 
	 * Constructs a continuous vector with points specified by a provided continuous vector
	 * 
	 * @param aCC	ContinuousVector which to initialise the points from
	 */
	public ContinuousVector(ContinuousVector aCC) {
		x = aCC.x;
		y = aCC.y;
		z = aCC.z;
	}

	/**
	 * \brief Constructs a continuous vector with points specified from XML tags
	 * 
	 * Constructs a continuous vector with points specified from XML tags
	 * 
	 * @param xmlRoot	Set of XML tags containing an X,Y,and Z coordinate
	 */
	public ContinuousVector(Element xmlRoot) {
		x = Double.parseDouble(xmlRoot.getAttributeValue("x"));
		y = Double.parseDouble(xmlRoot.getAttributeValue("y"));
		z = Double.parseDouble(xmlRoot.getAttributeValue("z"));
	}

	/**
	 * \brief Translate a discrete coordinates expressed on a discrete spatial grid with the resolution res to form continuous vector
	 * 
	 * Translate a discrete coordinates expressed on a discrete spatial grid with the resolution res to form continuous vector
	 * 
	 * @param dC	Discrete vector containing points on a grid
	 * @param res	The resolution of this grid, to use to transform these points
	 */
	public ContinuousVector(DiscreteVector dC, double res) {
		x = (.5+dC.i)*res;
		y = (.5+dC.j)*res;
		z = (.5+dC.k)*res;
	}

	/**
	 * \brief Create a continuous vector from three provided points
	 * 
	 * Create a continuous vector from three provided points
	 * 
	 * @param x	X coordinate
	 * @param y	Y coordinate
	 * @param z	Z coordinate
	 */
	public ContinuousVector(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/**
	 * \brief Set this vector to the points contained in a supplied continuous vector
	 * 
	 * Set this vector to the points contained in a supplied continuous vector
	 * 
	 * @param cc	Continuous vector of points to set this vector to
	 */
	public void set(ContinuousVector cc) {
		this.x = cc.x;
		this.y = cc.y;
		this.z = cc.z;
	}

	/**
	 * \brief Set this vector to the supplied X,Y,Z points
	 * 
	 * Set this vector to the supplied X,Y,Z points
	 * 
	 * @param x	X coordinate
	 * @param y Y coordinate
	 * @param z	Z coordinate
	 */
	public void set(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	
	/**
	 * \brief Set all points in the vector to zero
	 * 
	 * Set all points in the vector to zero
	 * 
	 */
	public void reset() {
		x = 0;
		y = 0;
		z = 0;
	}

	/**
	 * \brief Determine whether the numeric points in this vector are valid
	 * 
	 * Determine whether the numeric points in this vector are valid
	 * 
	 * @return	Boolean stating whether this vector is valid (true) or not (false)
	 */
	public boolean isValid() {
		return !(Double.isInfinite(x)|Double.isNaN(x)|Double.isInfinite(y)|Double.isNaN(y)
		        |Double.isInfinite(z)|Double.isNaN(z));
	}

	/**
	 * \brief Changes the sign of the vector. Used for movement vectors.
	 * 
	 * Changes the sign of the vector. Used for movement vectors.
	 */
	public void turnAround() {
		x = -x;
		y = -y;
		z = -z;
	}

	/**
	 * \brief Determine if this vector is in the location giving by the points X,Y,Z
	 * 
	 * Determine if this vector is in the location giving by the points X,Y,Z
	 * 
	 * @param x	X coordinate
	 * @param y Y coordinate
	 * @param z	Z coordinate
	 * @return Booloean stating whether the vector position and coordinate (x,y,z) is identical
	 */
	public boolean equals(double x, double y, double z) {
		if (this.x==x&this.y==y&this.z==z) return true;
		else return false;
	}

	/**
	 * \brief Determine if all points in the vector are zero
	 * 
	 * Determine if all points in the vector are zero
	 * 
	 * @return	Boolean stating whether all points in the vector are zero
	 */
	public boolean isZero() {
		return (this.x==0&this.y==0&&this.z==0);
	}

	/**
	 * \brief Print coordinates to string
	 * 
	 * Print coordinates to string
	 * 
	 * @return String containing the points in this vector
	 */
	public String toString() {
		return ExtraMath.toString(x, false)+",\t"+ExtraMath.toString(y, false)+",\t"
		        +ExtraMath.toString(z, false);
	}

	/**
	 * \brief Add vector v to this continuous vector
	 * 
	 * Add vector v to this continuous vector
	 * 
	 * @param v	ContinuousVector to add to this vector
	 */
	public void add(ContinuousVector v) {
		this.x += v.x;
		this.y += v.y;
		this.z += v.z;
	}

	/**
	 * \brief Add points X,Y,Z to their respective point in this vector
	 * 
	 * Add points X,Y,Z to their respective point in this vector
	 * 
	 * @param x	X coordinate
	 * @param y Y coordinate
	 * @param z	Z coordinater
	 */
	public void add(double x, double y, double z) {
		this.x += x;
		this.y += y;
		this.z += z;
	}

	/**
	 * \brief Store in this vector the sum of two other continuous vectors
	 * 
	 * Store in this vector the sum of two other continuous vectors
	 * 
	 * @param a	First continuous vector
	 * @param b	Continuous vector to add to first
	 */
	public void sendSum(ContinuousVector a, ContinuousVector b) 
	{
		this.x = a.x+b.x;
		this.y = a.y+b.y;
		this.z = a.z+b.z;
	}

	/**
	 * \brief Subtract vector v from this continuous vector
	 * 
	 * Subtract vector v from this continuous vector
	 * 
	 * @param v	ContinuousVector to subtract from this vector
	 */
	public void subtract(ContinuousVector v) {
		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;
	}

	/**
	 * \brief Store in this vector the subtraction of two other continuous vectors
	 * 
	 * Store in this vector the subtraction of two other continuous vectors
	 * 
	 * @param a	First continuous vector
	 * @param b	Continuous vector to subtract from the first
	 */
	public void sendDiff(ContinuousVector a, ContinuousVector b) {
		this.x = a.x-b.x;
		this.y = a.y-b.y;
		this.z = a.z-b.z;
	}

	/**
	 * \brief Calculate scalar product (dot product) of this vector with vector cc supplied
	 * 
	 * Calculate scalar product (dot product) of this vector with vector cc supplied
	 * 
	 * @param cc	Continuous vector to multiply (dot product) with this vector
	 * @return Double value of scalar product of two vectors
	 */
	public double prodScalar(ContinuousVector cc) {
		return this.x*cc.x+this.y*cc.y+this.z*cc.z;
	}

	/**
	 * \brief Multiply (stretch) this vector by supplied multiplier
	 * 
	 * Multiply (stretch) this vector by supplied multiplier
	 * 
	 * @param multiplier	Amount to stretch this vector by
	 */
	public void times(double multiplier) {
		this.x *= multiplier;
		this.y *= multiplier;
		this.z *= multiplier;
	}

	/**
	 * \brief Set this vector to a random coordinate (x,y,z)
	 * 
	 * Set this vector to a random coordinate (x,y,z)
	 * 
	 * @param is3D	Boolean noting if a Z coordinate needs to be calculated (if domain is 3D)
	 */
	public void alea(boolean is3D) {
		this.x = ExtraMath.getUniRand(-1,1);
		this.y = ExtraMath.getUniRand(-1,1);
		this.z = (is3D ? ExtraMath.getUniRand(-1,1) : 0);
	}

	/**
	 * \brief Normalize this Vector to unit length.
	 * 
	 * Normalize this Vector to unit length.
	 *
	 */
	public void normalizeVector() {
		double v = this.norm();
		if (v!=0) this.times(1/this.norm());

	}

	/**
	 * \brief Normalize this Vector to a given length.
	 * 
	 * Normalize this Vector to a given length.
	 * 
	 * @param newLength	Length used to normalise vector
	 */
	public void normalizeVector(double newLength) {
		this.times(newLength/this.norm());
	}

	/**
	 * \brief Calculate and return the absolute distance to a vector expressed in cc
	 * 
	 * Calculate and return the absolute distance to a vector expressed in cc
	 * 
	 * @param cc	ContinuousVector to calculate distance to
	 */
	public double distance(ContinuousVector cc) {
		return (double) Math.sqrt(Math.abs((this.x-cc.x))*Math.abs((this.x-cc.x))+Math.abs((this.y-cc.y))*Math.abs((this.y-cc.y))
		        + Math.abs((this.z-cc.z))*Math.abs((this.z-cc.z)));
	}

	/**
	 * \brief Return absolute length
	 * 
	 * Return absolute length
	 * 
	 * @return	Double value stating absolute length of this vector
	 */
	public double norm() {
		return Math.sqrt(x*x+y*y+z*z);
	}

	/**
	 * \brief Calculate cosine of the angle to vector cc
	 * 
	 * Calculate cosine of the angle to vector cc
	 * 
	 * @param cc	ContinuousVector for which cosine of the angle to this one should be calculated 	
	 * @return	Cosine of the angle to vector cc
	 */
	public double cosAngle(ContinuousVector cc) {
		return (x*cc.x+y*cc.y+z*cc.z)/Math.sqrt((x*x+y*y+z*z)*(cc.x*cc.x+cc.y*cc.y+cc.z*cc.z));
	}

	/**
	 * \brief Clone this vector, if supported
	 * 
	 * Clone this vector, if supported
	 * 
	 * @throws CloneNotSupportedException	Thrown if the object cannot be cloned
	 * 
	 */
	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
	
	/**
	 * \brief Compares the values of two continuous vectors. Used in the creation of the epithelium for eGUT
	 * 
	 * @param other	ContinuousVector to compare this continuous vector object to
	 */
	public int compareTo(ContinuousVector other) 
    {
		int valueComparison = Double.valueOf(this.x).compareTo(Double.valueOf(other.x));
    	
    	if (valueComparison == 0)
    		valueComparison = Double.valueOf(this.y).compareTo(Double.valueOf(other.y));
    	
    	if (valueComparison == 0)
    		valueComparison = Double.valueOf(this.z).compareTo(Double.valueOf(other.z));
    	
    	return valueComparison;
    }

}