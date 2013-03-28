/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ________________________________________
 * Abstract class with some extra math functions
 */

/**
 * @since May 2003
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */


package utils;

import java.text.DecimalFormat;
import java.util.*;

import Jama.Matrix;
import exceptions.ModelRuntimeException;

import simulator.geometry.ContinuousVector;

public final class ExtraMath {

	public static java.text.DecimalFormat dfSc = new DecimalFormat("000.###E0");
	public static java.text.DecimalFormat dfUs = new DecimalFormat("########.##");
	public static Random random;
	/**
	 * Computes the logarithm of base 2
	 * 
	 * @param x a number greater than 0.0
	 * @return the logarithm base 2 of x
	 */
	public static final double log2(double x) {
		return (Math.log(x)/Math.log(2));
	}

	/**
	 * Square
	 * 
	 * @param x value to square
	 * @return x*x
	 */
	public static final int sq(int x) {
		return x*x;
	}

	/**
	 * Square
	 * @param x: value to square
	 * @return x*x
	 */
	public static final double sq(double x) {
		return x*x;
	}

	/**
	 * cube of a number
	 * @param x : value to cube
	 * @return x*x*x
	 */
	public static final double cube(double x) {
		return x*x*x;
	}

	/**
	 * power of 2
	 * 
	 * @param x
	 * @return 2^x
	 */
	public static final int exp2(int x) {
		return (int) Math.pow(2, x);
	}

	/**
	 * power of 2
	 * 
	 * @param x
	 * @return 2^x
	 */
	public static final long exp2long(int x) {
		return (long) Math.pow(2, x);
	}
	
	/**
	 * power of 2
	 * 
	 * @param x
	 * @return 2^x
	 */
	public static final double exp2(double x) {
		return Math.pow(2, x);
	}
	
	/**
	 * The volume of a sphere with radius r
	 * @param r radius
	 * @return volume of sphere
	 */
	public static final double volumeOfASphere(double r) {
		return 4.1888d*r*r*r;
	}

	/**
	 * The volume of a cylinder with radius r and length l
	 * 
	 * @param r radius
	 * @param l length
	 * @return volume of cylinder
	 */
	public static final double volumeOfACylinder(double r, double l) {
		return 3.1416d*r*r*l;
	}

	/**
	 * The area of circle with radius r
	 * 
	 * @param r radius
	 * @return area of circle
	 */
	public static final double areaOfACircle(double r) {
		return 3.1416d*r*r;
	}

	/**
	 * Radius of a sphere with volume v
	 * 
	 * @param v volume
	 * @return radius
	 */
	public static final double radiusOfASphere(double v) {
		return (double) Math.pow(0.23873d*v, 0.33333d);
	}

	/**
	 * Radius of a sphere with volume v assuming repetition along 3rd dimension
	 * @param v volume
	 * @param Lz : omitted dimension
	 * @return radius
	 */
	public static final double radiusOfASphere(double v, double lZ) {
		return (double) Math.pow(0.23873d*v*lZ, 0.25d);
	}

	/**
	 * Returns the radius of a cylinder with volume v and length l
	 * 
	 * @param v
	 * @return
	 */
	public static final double radiusOfACylinder(double v, double l) {
		return (double) Math.sqrt(v/(3.1416d*l));
	}

	/**
	 * Distance between 2 points
	 * 
	 * @param x1,y1,z1 coordinates of point 1
	 * @param x2,y2,z2 coordinates of point 2
	 * @return distance
	 */
	public static final double pointDistance(double x1, double y1, double z1, double x2, double y2,
	        double z2) {
		return (double) Math.sqrt(sq(x1-x2)+sq(y1-y2)+sq(z1-z2));
	}

	/**
	 * Distance between 2 points
	 * 
	 * @param p1 point 1
	 * @param p2 point 2
	 * @return distance
	 */
	public static final double pointDistance(ContinuousVector p1, ContinuousVector p2) {
		return pointDistance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
	}

	/**
	 * Perform gamma correction of value v: v^gamma check that v is in the [0,
	 * 1] range
	 * 
	 * @param v
	 * @param gamma
	 * @return v^gamma
	 */
	public static double gammaCorrection(double v, double gamma) {
		if ((v<0)|(v>1)) throw new ModelRuntimeException("invalid v for gamma correction, v = "+v);
		return (double) Math.pow(v, gamma);
	}

	/**
	 * The maximum among 2 doubles
	 * 
	 * @param a
	 * @param b
	 * @return the maximum among a and b
	 */
	public static double max(double a, double b) {
		return (a>b ? a : b);
	}

	public static double min(double a, double b) {
		return (a<b ? a : b);
	}

	/**
	 * @param a Vector of double
	 * @return max double in a
	 */
	public static double max(double[] a) {
		double out = a[0];
		for (int i = 0; i<a.length; i++) {
			out = max(out, a[i]);
		}
		return out;
	}

	public static double average(double[] a) {
		double out=0;
		int n=0;
		for (int i = 0; i<a.length; i++)
			if(!Double.isInfinite(a[i])&&!Double.isNaN(a[i])){
					out += a[i];
					n++;
			}		
		return (out/n);
	}

	
	public static double min(double[] a) {
		double out = a[0];
		for (int i = 0; i<a.length; i++) {
			out = Math.min(out, a[i]);
		}
		return out;
	}

	/**
	 * The maximum square among 2 doubles
	 * 
	 * @param a
	 * @param b
	 * @return the maximum square among a and b
	 */
	public static double maxSquare(double a, double b) {
		double a2 = a*a;
		double b2 = b*b;
		return (a2>b2 ? a2 : b2);
	}

	static public boolean sameSign(double a, double b) {
		return (a*b>=0);
	}

	public static double sumVector(double vector[]) {
		double out = 0;
		for (int i = 0; i<vector.length; i++) {
			out += vector[i];
		}
		return out;
	}

	public static int moveX(int index, int[] gridDim, int sens) {
		return index+sens;
	}

	public static int moveY(int index, int[] gridDim, int sens) {
		return index+sens*gridDim[0];
	}

	public static int moveZ(int index, int[] gridDim, int sens) {
		return index+sens*gridDim[0]*gridDim[1];
	}
	
	public static double truncateDigit(double value, int digits) {
		double n = Math.floor(Math.log10(Math.abs(value)));
		value = Math.round(value/Math.pow(10, n-digits-1));
		return value*Math.pow(10, n-digits-1);
	}

	public static String toString(double value, boolean scFormat) {
		if (scFormat) return dfSc.format(value);
		else return dfUs.format(value);
	}
	
	/**
	 * @return uniform distributed random number in [0,1]
	 */
	public static double getUniRand() {
		return random.nextDouble();
	}

	/**
	 * @param lBound
	 * @param hBound
	 * @return an int between lBound inclusive and hBound exclusive
	 */
	public static int getUniRandInt(int lBound, int hBound) {
		return random.nextInt(hBound-lBound)+lBound;
	}
	
	/**
	 * @param lBound
	 * @param hBound
	 * @return a double between lBound inclusive and hBound exclusive
	 */
	public static double getUniRand(double lBound, double hBound) {
		return random.nextDouble()*(hBound-lBound)+lBound;
	}

	/**
	 * @return truncated N(0,1) distributed random number. Normal distributed
	 * random numbers are truncated at 2*sigma to prevent extreme values.
	 */
	public static double getNormRand() {
		double phi;
		do {
			phi = random.nextGaussian();
		} while (Math.abs(phi)>2);
		return phi;
	}
	
	/**
	 * @return 2 to the power of a uniformly distributed random number in [0,1]
	 */
	public static double getExp2Rand() {
		return exp2(getUniRand());
	}
	
	/**
	 * Randomise a value with a gaussian distribution in a range fixed by the
	 * ratio sigma. If mu = 1 and sigma = .1, the results form a truncated
	 * gaussian distribution between 0.8 and 1.2 (2*sigma interval)
	 * 
	 * @param mu mean value
	 * @param sigma standard deviation
	 * @return N(mu,sigma)-distributed random value within [mu-2*sigma,
	 * mu+2*sigma]
	 */
	public static double deviateFrom(double mu, double sigma) {
		// Called by getBabyMassFrac(), willDie(), willDivide() and willTransfer().
		// sigma is homogeneous to a ratio of deviation/value
		
		// Rob 11th Jan 2011 
		// Extended "if (mu==0) return 0;" to avoid unnecessarily calling getNormRand()
		if (mu==0 || sigma==0) return mu;
		
		double result, deviation;

		do {
			deviation = getNormRand();
			// Compute the new value
			result = mu*(1.0+sigma*deviation);
		} while (!sameSign(result, mu));

		return result;
	}
	
	// bvm 16.12.08
	public static double mean(double [] v) {
		return sumVector(v)/v.length;
	}

	// bvm 16.12.08
	public static double stddev(double [] v) {
		double mean = mean(v);
		double sum = 0.;
		for (int i = 0; i<v.length; i++) {
			sum += sq(v[i] - mean);
		}

		return Math.sqrt(sum/v.length);
	}

	                                      
}
