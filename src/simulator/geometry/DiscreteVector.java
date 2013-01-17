/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *______________________________________________________
 * Implements 3D vector of discrete spatial coordinates
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */


package simulator.geometry;

import java.io.Serializable;

import org.jdom.Element;


public class DiscreteVector implements Cloneable, Serializable {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	public int i;
	public int j;
	public int k;

	/**
     * Constructor (empty vector)
     */
	public DiscreteVector() {
		i = 0;
		j = 0;
		k = 0;
	}
	
	/**
	 * XML-based constructor
	 * @param coordinatesRoot
	 */
	public DiscreteVector(Element coordinatesRoot){
		i=Integer.parseInt(coordinatesRoot.getAttributeValue("x"));
		j=Integer.parseInt(coordinatesRoot.getAttributeValue("y"));
		k=Integer.parseInt(coordinatesRoot.getAttributeValue("z"));
	}

	public DiscreteVector clone() {
		DiscreteVector out = new DiscreteVector(i, j, k);
		return out;
	}

	public void set(DiscreteVector dV) {
		this.i = dV.i;
		this.j = dV.j;
		this.k = dV.k;
	}
	
	/**
     * Translate a continuous coordinates to a discrete one on a discrete
     * spatial grid with the resolution res
     * @param cc
     * @param res grid resolution
     */
	public DiscreteVector(ContinuousVector cc, double res) {
		i = (int) Math.ceil(cc.x/res);
		j = (int) Math.ceil(cc.y/res);
		k = (int) Math.ceil(cc.z/res);
	}

	public void set(int i0, int j0, int k0) {
		this.i = i0;
		this.j = j0;
		this.k = k0;

	}

	/**
     * set coordinate to origin
     */
	public void reset() {
		i = 0;
		j = 0;
		k = 0;
	}

	/**
     * Constructor
     * @param l depth coordinate
     * @param m horizontal coordinate
     * @param n vertical coordinate
     */
	public DiscreteVector(int n, int m, int l) {
		this.k = l;
		this.j = m;
		this.i = n;
	}

	public void add(int i, int j, int k) {
		this.i += i;
		this.j += j;
		this.k += k;
	}

	public void add(DiscreteVector dC) {
		this.i += dC.i;
		this.j += dC.j;
		this.k += dC.k;
	}

	public void sendSum(DiscreteVector a, DiscreteVector b) {
		i = a.i+b.i;
		j = a.j+b.j;
		k = a.k+b.k;
	}

	public void diff(DiscreteVector dC) {
		this.i -= dC.i;
		this.j -= dC.j;
		this.k -= dC.k;
	}

	public void times(double n) {
		this.i *= n;
		this.j *= n;
		this.k *= n;
	}

	public void turnAround() {
		i = -i;
		j = -j;
		k = -k;
	}

	public double norm() {
		return Math.sqrt(i*i+j*j+k*k);
	}

	public boolean equals(DiscreteVector dc) {
		return ((dc.i==this.i)&(dc.j==this.j)&(dc.k==this.k));
	}

	public int prodScalar(DiscreteVector dc) {
		return this.i*dc.i+this.j*dc.j+this.k*dc.k;
	}
	
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
	/*
     * public Object clone() throws CloneNotSupportedException { return
     * super.clone(); }
     */
	public String toString(){
		return "("+i+"-"+j+"-"+k+")";
	}
}
