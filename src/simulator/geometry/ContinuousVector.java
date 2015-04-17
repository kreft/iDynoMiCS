/**
 * \package simulator.geometry
 * \brief Package of boundary utilities that aid the creation of the
 * environment being simulated.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.geometry;

import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Implements 3D vector of continuous spatial coordinates.
 * 
 * Cartesian (x, y, z) coordinates obligatory.
 * Can be used to store Continuous coordinates or Movement vectors.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA).
 *
 */
public class ContinuousVector implements Cloneable
{
	/**
	 * X coordinate of the point contained in this vector
	 */
	public Double x;
	
	/**
	 * Y coordinate of the point contained in this vector
	 */
	public Double y; 
	
	/**
	 * Z coordinate of the point contained in this vector
	 */
	public Double z;
	
	/**
	 * \brief Constructs a ContinuousVector at origin.
	 */
	public ContinuousVector()
	{
		reset();
	}
	
	/**
	 * \brief Constructs a continuous vector with points specified by a
	 * provided continuous vector.
	 * 
	 * @param aCC ContinuousVector which to initialise the points from.
	 */
	public ContinuousVector(ContinuousVector aCC)
	{
		set(aCC);
	}
	
	/**
	 * \brief Constructs a continuous vector with points specified from an
	 * XMLParser.
	 * 
	 * @param coordinatesRoot	An XMLParser containing x, y and z coordinates.
	 */
	public ContinuousVector(XMLParser coordinatesRoot)
	{
		x = coordinatesRoot.getAttributeDbl("x");
		y = coordinatesRoot.getAttributeDbl("y");
		z = coordinatesRoot.getAttributeDbl("z");
	}
	
	/**
	 * \brief Create a continuous vector from three provided points.
	 * 
	 * @param x	X coordinate.
	 * @param y	Y coordinate.
	 * @param z	Z coordinate.
	 */
	public ContinuousVector(Double x, Double y, Double z)
	{
		set(x, y, z);
	}
	
	/**
	 * \brief Set this vector to the points contained in a supplied continuous
	 * vector.
	 * 
	 * @param cc Continuous vector of points to set this vector to.
	 */
	public void set(ContinuousVector cc)
	{
		set(cc.x, cc.y, cc.z);
	}
	
	/**
	 * \brief Set this vector to the supplied X,Y,Z points.
	 * 
	 * @param x	X coordinate.
	 * @param y Y coordinate.
	 * @param z	Z coordinate.
	 */
	public void set(Double x, Double y, Double z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	/**
	 * 
	 * TODO Rob 13Mar2015: Check padding.
	 * Compare Agentcontainer.getGridLocation(int index)
	 * 
	 * @param dC
	 * @param res
	 */
	public void setToVoxelCenter(DiscreteVector dC, Double res)
	{
		set(dC.i + 0.5, dC.j + 0.5, dC.k + 0.5);
		times(res);
	}
	
	/**
	 * 
	 * @param dC
	 */
	public void set(DiscreteVector dC)
	{
		set(dC.i + 0.0, dC.j + 0.0, dC.k + 0.0);
	}
	
	/**
	 * 
	 * @param dC
	 * @param res
	 */
	public void set(DiscreteVector dC, Double res)
	{
		set(dC);
		times(res);
	}
	
	/**
	 * \brief Set all points in the vector to zero.
	 */
	public void reset()
	{
		set(0.0, 0.0, 0.0);
	}

	/**
	 * \brief Determine whether the numeric points in this vector are valid.
	 * 
	 * @return Boolean stating whether this vector is valid (true) or not (false)
	 */
	public Boolean isValid()
	{
		return ! ( x.isInfinite() || x.isNaN() ||
				   y.isInfinite() || y.isNaN() ||
				   z.isInfinite() || z.isNaN() );
	}

	/**
	 * \brief Changes the sign of the vector.
	 * 
	 * Used for movement vectors.
	 */
	public void turnAround()
	{
		times(-1.0);
	}

	/**
	 * \brief Determine if this vector is in the given location.
	 * 
	 * @param x	X coordinate.
	 * @param y Y coordinate.
	 * @param z	Z coordinate.
	 * @return Boolean stating whether the vector position and coordinate
	 * (x,y,z) is identical.
	 */
	public Boolean equals(Double x, Double y, Double z)
	{
		return (this.x.equals(x) && this.y.equals(y) && this.z.equals(z));
	}
	
	/**
	 * \brief Check if this vector is the same as another.
	 *  
	 * @param other
	 * @return
	 */
	public Boolean equals(ContinuousVector other)
	{
		return equals(other.x, other.y, other.z);
	}
	
	/**
	 * \brief Determine if all points in the vector are zero.
	 * 
	 * @return	Boolean stating whether all points in the vector are zero.
	 */
	public Boolean isZero()
	{
		return equals(0.0, 0.0, 0.0);
	}
	
	/**
	 * \brief Print coordinates to string.
	 * 
	 * @return String containing the points in this vector.
	 */
	@Override
	public String toString()
	{
		return "("+ExtraMath.toString(x, false)
				+", "+ExtraMath.toString(y, false)
				+", "+ExtraMath.toString(z, false)+")";
	}
	
	/**
	 * \brief Add vector v to this continuous vector.
	 * 
	 * @param v	ContinuousVector to add to this vector.
	 */
	public void add(ContinuousVector v)
	{
		add(v.x, v.y, v.z);
	}
	
	/**
	 * \brief Add points X,Y,Z to their respective point in this vector.
	 * 
	 * @param x	X coordinate.
	 * @param y Y coordinate.
	 * @param z	Z coordinate.
	 */
	public void add(Double x, Double y, Double z)
	{
		this.x += x;
		this.y += y;
		this.z += z;
	}
	
	/**
	 * \brief Store in this vector the sum of two other continuous vectors.
	 * 
	 * @param a	First continuous vector.
	 * @param b	Continuous vector to add to first.
	 */
	public void sendSum(ContinuousVector a, ContinuousVector b) 
	{
		set(a);
		add(b);
	}
	
	/**
	 * \brief Subtract vector v from this continuous vector.
	 * 
	 * @param cV	ContinuousVector to subtract from this vector.
	 */
	public void subtract(ContinuousVector cV)
	{
		add(-cV.x, -cV.y, -cV.z);
	}

	/**
	 * \brief Store in this vector the difference of two other continuous
	 * vectors.
	 * 
	 * @param a	First continuous vector.
	 * @param b	Continuous vector to subtract from the first.
	 */
	public void sendDiff(ContinuousVector a, ContinuousVector b)
	{
		set(a);
		subtract(b);
	}
	
	/**
	 * \brief Calculate scalar product (dot product) of this vector with vector
	 * cc supplied.
	 * 
	 * @param cc Continuous vector to multiply (dot product) with this vector.
	 * @return Double value of scalar product of two vectors.
	 */
	public Double prodScalar(ContinuousVector cc)
	{
		return (this.x * cc.x) + (this.y * cc.y) + (this.z * cc.z);
	}

	/**
	 * \brief Multiply (stretch) this vector by supplied multiplier.
	 * 
	 * @param multiplier Amount to stretch this vector by.
	 */
	public void times(Double multiplier)
	{
		this.x *= multiplier;
		this.y *= multiplier;
		this.z *= multiplier;
	}

	/**
	 * \brief Set this vector to a random coordinate (x,y,z). 
	 * 
	 * Coordinates are chosen from a uniform distribution in (-1, 1). 
	 * 
	 * @param is3D Boolean noting if a Z coordinate needs to be calculated.
	 */
	public void alea(Boolean is3D)
	{
		this.x = ExtraMath.getUniRand(-1.0, 1.0);
		this.y = ExtraMath.getUniRand(-1.0, 1.0);
		this.z = (is3D ? ExtraMath.getUniRand(-1.0, 1.0) : 0.0);
	}
	
	/**
	 * \brief Normalize this Vector to unit length.
	 */
	public void normalizeVector()
	{
		normalizeVector(1.0);
	}

	/**
	 * \brief Normalize this Vector to a given length.
	 * 
	 * @param newLength	Length used to normalise vector.
	 */
	public void normalizeVector(Double newLength)
	{
		Double norm = this.norm();
		if ( ! norm.equals(0.0) )
			this.times(newLength/norm);
	}

	/**
	 * \brief Calculate and return the absolute distance to a vector expressed in cc.
	 * 
	 * Does not take cyclic boundaries into account. 
	 * 
	 * @param cc ContinuousVector to calculate distance to.
	 */
	public Double distance(ContinuousVector cc)
	{
		return ExtraMath.hypotenuse(this.x - cc.x, this.y - cc.y, this.z - cc.z);
	}

	/**
	 * \brief Return absolute length.
	 * 
	 * @return	Double value stating absolute length of this vector.
	 */
	public Double norm()
	{
		return ExtraMath.hypotenuse(x, y, z);
	}
	
	/**
	 * \brief Calculate cosine of the angle to a given vector.
	 * 
	 * @param v ContinuousVector for which cosine of the angle to this one
	 * should be calculated.
	 * @return	Cosine of the angle to vector given.
	 */
	public Double cosAngle(ContinuousVector v)
	{
		/*
		 * Returning 0.0 if the dot product is 0.0 removes the danger of
		 * dividing 0.0/0.0 and also speeds things up slightly if the vectors
		 * are orthogonal.
		 */
		Double dotProd = prodScalar(v);
		return ( dotProd == 0.0 ) ? 0.0 : dotProd/(this.norm() * v.norm());
	}
	
	/**
	 * \brief Calculate the angle to another vector.
	 * 
	 * @param v	Another vector
	 * @return Angle between this vector and that given.
	 */
	public Double angle(ContinuousVector v)
	{
		return Math.acos(cosAngle(v));
	}
	
	/**
	 * \brief Clone this vector, if supported.
	 * 
	 * @throws CloneNotSupportedException Thrown if the object cannot be cloned.
	 */
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		return super.clone();
	}
	
	public Boolean isOrthogonal(ContinuousVector other)
	{
		return ( prodScalar(other) == 0.0 );
	}
	
	public Boolean isParallel(ContinuousVector other)
	{
		return ( Math.abs(cosAngle(other)) == 1.0 );
	}
	
	
	/**
	 * \brief Returns the cross product of this vector with another.
	 * 
	 * Note that the resulting vector is orthogonal to both this vector and
	 * the one given.
	 * 
	 * @param other
	 * @return
	 */
	public ContinuousVector crossProduct(ContinuousVector other)
	{
		ContinuousVector out = new ContinuousVector();
		out.x = this.y*other.z - this.z*other.y;
		out.y = this.z*other.x - this.x*other.z;
		out.z = this.x*other.y - this.y*other.x;
		return out;
	}
}