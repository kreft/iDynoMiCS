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

import org.jdom.Element;

import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Implements 3D vector of discrete spatial coordinates.
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA).
 */
public class DiscreteVector implements Cloneable, Serializable 
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * I Location on a grid.
	 */
	public int i;
	
	/**
	 * J Location on a grid.
	 */
	public int j;
	
	/**
	 * K Location on a grid.
	 */
	public int k;
	
	/**
     * \brief Creates a discrete vector initialised at the origin.
     */
	public DiscreteVector()
	{
		reset();
	}
	
	/**
	 * 
	 * @param dV
	 */
	public DiscreteVector(DiscreteVector dV)
	{
		set(dV);
	}
	
	/**
	 * \brief Constructs a discrete vector with points specified from an
	 * XMLParser.
	 * 
	 * @param coordinatesRoot	An XMLParser containing i, j and k coordinates.
	 */
	public DiscreteVector(XMLParser coordinatesRoot)
	{
		i = coordinatesRoot.getAttributeInt("i");
		j = coordinatesRoot.getAttributeInt("j");
		k = coordinatesRoot.getAttributeInt("k");
		
		try
		{
			i = coordinatesRoot.getAttributeInt("i");
			j = coordinatesRoot.getAttributeInt("j");
			k = coordinatesRoot.getAttributeInt("k");
		}
		catch (Exception e)
		{
			i = (int) (double) coordinatesRoot.getAttributeDbl("x");
			j = (int) (double) coordinatesRoot.getAttributeDbl("y");
			k = (int) (double) coordinatesRoot.getAttributeDbl("z");
			LogFile.writeLogAlways("----------------------------------------");
			LogFile.writeLogAlways("Please give discrete vectors as i, j, k!");
			LogFile.writeLogAlways("x, y, z should be for continuous vectors");
			LogFile.writeLogAlways("----------------------------------------");
		}
		
	}
	
	/**
	 * \brief Creates a clone of this discrete vector.
	 * 
	 * @return Clone of this discrete vector object.
	 */
	@Override
	public DiscreteVector clone()
	{
		DiscreteVector out = new DiscreteVector(i, j, k);
		return out;
	}
	
	/**
	 * \brief Constructs a continuous vector with points specified by a
	 * provided discrete vector.
	 * 
	 * @param dV DiscreteVector from which to initialise the points.
	 */
	public void set(DiscreteVector dV)
	{
		set(dV.i, dV.j, dV.k);
	}
	
	/**
	 * \brief Translate a continuous coordinate expressed on a spatial grid
	 * with the resolution res to form a discrete vector.
	 * 
	 * @param cV Continuous vector containing points on a grid.
	 * @param res The resolution of this grid, to use to transform these points.
	 */
	public DiscreteVector(ContinuousVector cV, Double res)
	{
		set(cV, res);
	}
	
	/**
	 * \brief Translate a continuous coordinate expressed on a spatial grid
	 * with the resolution res to form a discrete vector.
	 * 
	 * @param cV Continuous vector containing points on a grid.
	 * @param res The resolution of this grid, to use to transform these points.
	 */
	public void set(ContinuousVector cV, Double res)
	{
		set( (int) Math.ceil(cV.x/res),
			 (int) Math.ceil(cV.y/res),
			 (int) Math.ceil(cV.z/res));
	}
	
	/**
	 * \brief Set this vector to the supplied i, j, k points.
	 * 
	 * @param i i-coordinate
	 * @param j j-coordinate
	 * @param k k-coordinate
	 */
	public void set(int i, int j, int k)
	{
		this.i = i;
		this.j = j;
		this.k = k;
	}

	/**
	 * \brief Set all points in the vector to zero.
	 */
	public void reset()
	{
		set(0, 0, 0);
	}

	/**
	 * \brief Create a discrete vector from the supplied i, j, k points.
	 * 
	 * @param i i-coordinate
	 * @param j j-coordinate
	 * @param k k-coordinate
	 */
	public DiscreteVector(int i, int j, int k)
	{
		set(i, j, k);
	}
	
	/**
	 * \brief Add points I,J,K to their respective point in this vector.
	 * 
	 * @param i	I coordinate.
	 * @param j J coordinate.
	 * @param k	K coordinate.
	 */
	public void add(int i, int j, int k)
	{
		this.i += i;
		this.j += j;
		this.k += k;
	}

	/**
	 * \brief Add vector v to this discrete vector.
	 * 
	 * @param dV DiscreteVector to add to this vector.
	 */
	public void add(DiscreteVector dV)
	{
		add(dV.i, dV.j, dV.k);
	}
	
	/**
	 * \brief Store in this vector the sum of two other discrete vectors.
	 * 
	 * @param a	First discrete vector.
	 * @param b	Discrete vector to add to first.
	 */
	public void sendSum(DiscreteVector a, DiscreteVector b)
	{
		set(a);
		add(b);
	}
	
	/**
	 * \brief Subtract vector v from this discrete vector.
	 * 
	 * @param dV DiscreteVector to subtract from this vector.
	 */
	public void subtract(DiscreteVector dV) 
	{
		add( -dV.i, -dV.j, -dV.k );
	}
	
	/**
	 * \brief Store in this vector the difference of two other discrete
	 * vectors.
	 * 
	 * @param a	First discrete vector.
	 * @param b	Discrete vector to subtract from the first.
	 */
	public void sendDiff(DiscreteVector a, DiscreteVector b)
	{
		set(a);
		subtract(b);
	}
	
	/**
	 * \brief Multiply (stretch) this vector by supplied multiplier.
	 * 
	 * @param n	Amount to stretch this vector by.
	 */
	public void times(int n)
	{
		i *= n;
		j *= n;
		k *= n;
	}

	/**
	 * \brief Changes the sign of the vector.
	 * 
	 * Used for movement vectors.
	 */
	public void turnAround()
	{
		times(-1);
	}

	/**
	 * \brief Return absolute length of this vector.
	 * 
	 * @return Double value stating absolute length of this vector.
	 */
	public Double norm()
	{
		return ExtraMath.hypotenuse(i, j, k);
	}
	
	/**
	 * \brief Calculate cosine of the angle to a given vector.
	 * 
	 * @param v DiscreteVector for which cosine of the angle to this one
	 * should be calculated.
	 * @return	Cosine of the angle to vector given.
	 */
	public Double cosAngle(DiscreteVector v)
	{
		return prodScalar(v)/(this.norm() * v.norm());
	}
	
	public Boolean isParallel(DiscreteVector other)
	{
		return ( Math.abs(cosAngle(other)) == 1.0 );
	}
	
	/**
	 * \brief Determine if this vector equals the points given in the provided
	 * vector.
	 * 
	 * @param dc Discrete vector to compare to this vector
	 * @return	Boolean stating whether the two vectors are equal
	 */
	public Boolean equals(DiscreteVector dc)
	{
		return equals(dc.i, dc.j, dc.k);
	}
	
	/**
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	public Boolean equals(int i, int j, int k)
	{
		return ( (i == this.i) && (j == this.j) && (k == this.k));
	}
	
	/**
	 * 
	 * @return
	 */
	public Boolean isZero()
	{
		return equals(0, 0, 0);
	}
	
	/**
	 * \brief Calculate scalar product (dot product) of this vector with the
	 * discrete vector dV supplied.
	 * 
	 * @param dV The discrete vector to multiply (dot product) with this vector.
	 * @return Integer value of scalar product of two vectors.
	 */
	public int prodScalar(DiscreteVector dV)
	{
		return (this.i * dV.i) + (this.j * dV.j) + (this.k * dV.k);
	}
	
	/**
	 * \brief Calculates two orthogonal vectors colinear to this vector.
	 * 
	 * @param v	First discrete vector to produce.
	 * @param w Second discrete vector to produce.
	 */
	public void orthoVector(DiscreteVector v, DiscreteVector w)
	{
		if ( i != 0 )
		{
			v.set(-j/i, 1, 0);
			if ( k != 0 )
				w.set(1, j/i, -( i + j*j/i )/k);
			else 
				w.set(0, 0, 1);
		}
		else if ( j != 0 )
		{
			v.set(0, -k/j, 1);
			w.set(1, 0, 0);
		}
		else
		{
			v.set(1, 0, 0);
			w.set(0, 1, 0);
		}
	}
	
	/**
	 * \brief Print coordinates to string.
	 * 
	 * @return String containing the points in this vector.
	 */
	@Override
	public String toString()
	{
		return "("+i+"-"+j+"-"+k+")";
	}
}
