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

/**
 * \brief Implements 3D vector of discrete spatial coordinates
 * 
 * Implements 3D vector of discrete spatial coordinates
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public class DiscreteVector implements Cloneable, Serializable 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * I Location on a grid
	 */
	public int i;
	
	/**
	 * J Location on a grid
	 */
	public int j;
	
	/**
	 * K Location on a grid
	 */
	public int k;

	/**
     * \brief Creates a discrete vector initialised at 0,0,0
     * 
     * Creates a discrete vector initialised at 0,0,0
     */
	public DiscreteVector() {
		i = 0;
		j = 0;
		k = 0;
	}
	
	/**
	 * \brief Constructs a discrete vector with points specified from XML tags
	 * 
	 * Constructs a discrete vector with points specified from XML tags
	 * 
	 * @param coordinatesRoot	Set of XML tags containing an X,Y,and Z coordinate
	 */
	public DiscreteVector(Element coordinatesRoot){
		i=Integer.parseInt(coordinatesRoot.getAttributeValue("x"));
		j=Integer.parseInt(coordinatesRoot.getAttributeValue("y"));
		k=Integer.parseInt(coordinatesRoot.getAttributeValue("z"));
	}

	/**
	 * \brief Creates a clone of this discrete vector
	 * 
	 * Creates a clone of this discrete vector
	 * 
	 * @return Clone of this discrete vector object
	 */
	public DiscreteVector clone() {
		DiscreteVector out = new DiscreteVector(i, j, k);
		return out;
	}

	/**
	 * \brief Constructs a continuous vector with points specified by a provided discrete vector
	 * 
	 * Constructs a continuous vector with points specified by a provided discrete vector
	 * 
	 * @param dV	DiscreteVector which to initialise the points from
	 */
	public void set(DiscreteVector dV) {
		this.i = dV.i;
		this.j = dV.j;
		this.k = dV.k;
	}
	
	/**
	 * \brief Translate a continuous coordinate expressed on a spatial grid with the resolution res to form a discrete vector
	 * 
	 * Translate a continuous coordinate expressed on a spatial grid with the resolution res to form a discrete vector
	 * 
	 * @param cc	Continuous vector containing points on a grid
	 * @param res	The resolution of this grid, to use to transform these points
	 */
	public DiscreteVector(ContinuousVector cc, double res) {
		i = (int) Math.ceil(cc.x/res);
		j = (int) Math.ceil(cc.y/res);
		k = (int) Math.ceil(cc.z/res);
	}

	/**
	 * \brief Set this vector to the supplied i,j,k points
	 * 
	 * Set this vector to the supplied i,j,k points
	 * 
	 * @param i0	i coordinate
	 * @param j0 j coordinate
	 * @param k0	k coordinate
	 */
	public void set(int i0, int j0, int k0) {
		this.i = i0;
		this.j = j0;
		this.k = k0;

	}

	/**
	 * \brief Set all points in the vector to zero
	 * 
	 * Set all points in the vector to zero
	 * 
	 */
	public void reset() {
		i = 0;
		j = 0;
		k = 0;
	}

	/**
	 * \brief Create a discrete vector from three provided points
	 * 
	 * Create a discrete vector from three provided points
	 * 
	 * @param n	N coordinate
	 * @param m	M coordinate
	 * @param l	L coordinate
	 */
	public DiscreteVector(int n, int m, int l) {
		this.k = l;
		this.j = m;
		this.i = n;
	}

	/**
	 * \brief Add points I,J,K to their respective point in this vector
	 * 
	 * Add points I,J,K to their respective point in this vector
	 * 
	 * @param i	I coordinate
	 * @param j J coordinate
	 * @param k	K coordinate
	 */
	public void add(int i, int j, int k) {
		this.i += i;
		this.j += j;
		this.k += k;
	}

	/**
	 * \brief Add vector v to this discrete vector
	 * 
	 * Add vector v to this discrete vector
	 * 
	 * @param dC	DiscreteVector to add to this vector
	 */
	public void add(DiscreteVector dC) {
		this.i += dC.i;
		this.j += dC.j;
		this.k += dC.k;
	}

	/**
	 * \brief Store in this vector the sum of two other discrete vectors
	 * 
	 * Store in this vector the sum of two other discrete vectors
	 * 
	 * @param a	First discrete vector
	 * @param b	Discrete vector to add to first
	 */
	public void sendSum(DiscreteVector a, DiscreteVector b) {
		i = a.i+b.i;
		j = a.j+b.j;
		k = a.k+b.k;
	}

	/**
	 * \brief Subtract vector v from this discrete vector
	 * 
	 * Subtract vector v from this discrete vector
	 * 
	 * @param dC	DiscreteVector to subtract from this vector
	 */
	public void diff(DiscreteVector dC) 
	{
		this.i -= dC.i;
		this.j -= dC.j;
		this.k -= dC.k;
	}

	/**
	 * \brief Multiply (stretch) this vector by supplied multiplier
	 * 
	 * Multiply (stretch) this vector by supplied multiplier
	 * 
	 * @param n	Amount to stretch this vector by
	 * 
	 */
	public void times(double n) {
		this.i *= n;
		this.j *= n;
		this.k *= n;
	}

	/**
	 * \brief Changes the sign of the vector. Used for movement vectors.
	 * 
	 * Changes the sign of the vector. Used for movement vectors.
	 */
	public void turnAround() {
		i = -i;
		j = -j;
		k = -k;
	}

	/**
	 * \brief Return absolute length
	 * 
	 * Return absolute length
	 * 
	 * @return	Double value stating absolute length of this vector
	 */
	public double norm() {
		return Math.sqrt(i*i+j*j+k*k);
	}

	/**
	 * \brief Determine if this vector equals the points given in the provided vector
	 * 
	 * Determine if this vector equals the points given in the provided vector
	 * 
	 * @param dc	Discrete vector to compare to this vector
	 * @return	Boolean stating whether the two vectors are equal
	 */
	public boolean equals(DiscreteVector dc) {
		return ((dc.i==this.i)&(dc.j==this.j)&(dc.k==this.k));
	}

	/**
	 * \brief Calculate scalar product (dot product) of this vector with vector dc supplied
	 * 
	 * Calculate scalar product (dot product) of this vector with vector dc supplied
	 * 
	 * @param dc	Discrete vector to multiply (dot product) with this vector
	 * @return Double value of scalar product of two vectors
	 */
	public int prodScalar(DiscreteVector dc) {
		return this.i*dc.i+this.j*dc.j+this.k*dc.k;
	}
	
	/**
	 * \brief Calculates two orthogonal vectors colinear to this vector
	 * 
	 * Calculates two orthogonal vectors colinear to this vector
	 * 
	 * @param v	First discrete vector to produce
	 * @param w Second discrete vector to produce
	 */
	public void orthoVector(DiscreteVector v, DiscreteVector w) {
		if (this.i!=0) {
			v.i = -j/i;
			v.j = 1;
			v.k = 0;
			if (this.k!=0) {
				w.i = 1;
				w.j = j/i;
				w.k = -(i+j*j/i)/k;
			} else {
				w.i = 0;
				w.j = 0;
				w.k = 1;
			}
		} else if (this.j!=0) {
			v.i = 0;
			v.j = -k/j;
			v.k = 1;
			w.i = 1;
			w.j = 0;
			w.k = 0;
		} else {
			v.i = 1;
			v.j = 0;
			v.k = 0;
			
			w.i = 0;
			w.j = 1;
			w.k = 0;
		}

	}
	
	/**
	 * \brief Print coordinates to string
	 * 
	 * Print coordinates to string
	 * 
	 * @return String containing the points in this vector
	 */
	public String toString(){
		return "("+i+"-"+j+"-"+k+")";
	}
}
