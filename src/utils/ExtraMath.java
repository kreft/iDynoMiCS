/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import java.text.DecimalFormat;
import java.util.*;

import Jama.Matrix;
import exceptions.ModelRuntimeException;

import simulator.geometry.ContinuousVector;

/**
 * \brief Abstract class with some extra useful math functions
 * 
 * Abstract class with some extra useful math functions
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public final class ExtraMath 
{

	/**
	 * One option for decimal of decimals when these are written to the screen
	 */
	public static java.text.DecimalFormat dfSc = new DecimalFormat("000.###E0");
	
	/**
	 * Second option for decimal of decimals when these are written to the screen
	 */
	public static java.text.DecimalFormat dfUs = new DecimalFormat("########.##");
	
	/**
	 * Random number generator
	 */
	public static Random random;
	
	/**
	 * \brief Computes the logarithm of base 2
	 * 
	 * Computes the logarithm of base 2
	 * 
	 * @param x a number greater than 0.0
	 * @return the logarithm base 2 of x
	 */
	public static final double log2(double x) {
		return (Math.log(x)/Math.log(2));
	}

	/**
	 * \brief Square an integer number
	 * 
	 * Square an integer number
	 * 
	 * @param x value to square
	 * @return The square of x
	 */
	public static final int sq(int x) {
		return x*x;
	}

	/**
	 * \brief Square a double
	 * 
	 * Square a double
	 * 
	 * @param x: value to square
	 * @return The square value of x
	 */
	public static final double sq(double x) {
		return x*x;
	}

	/**
	 * \brief cube a double number
	 * 
	 * cube a double number
	 * 
	 * @param x : value to cube
	 * @return The cube of x
	 */
	public static final double cube(double x) {
		return x*x*x;
	}

	/**
	 * \brief Calculate x to the power of x where x is an integer, returning a long number
	 * 
	 * Calculate x to the power of x where x is an integer, returning a long number
	 * 
	 * @param x	Number to be raised to power of 2
	 * @return	2^x, as a long number
	 */
	public static final long exp2long(int x) 
	{
		return (long) Math.pow(2, x);
	} 
	
	/**
	 * \brief Calculate x to the power of 2 where x is an integer
	 * 
	 * Calculate x to the power of 2 where x is an integer
	 * 
	 * @param x Number to be raised to power of 2
	 * @return 2^x
	 */
	public static final int exp2(int x) {
		return (int) Math.pow(2, x);
	}

	/**
	 * \brief Calculate x to the power of 2 where x is a double
	 * 
	 * Calculate x to the power of 2 where x is an double
	 * 
	 * @param x Number to be raised to power of 2
	 * @return 2^x
	 */
	public static final double exp2(double x) {
		return Math.pow(2, x);
	}
	
	/**
	 * \brief Calculate the volume of a sphere with radius r
	 * 
	 * Calculate the volume of a sphere with radius r
	 * 
	 * @param r	Radius of the sphere
	 * @return volume of sphere
	 */
	public static final double volumeOfASphere(double r) {
		return 4.1888d*r*r*r;
	}

	/**
	 * \brief Calculate the volume of a cylinder with radius r and length l
	 * 
	 * Calculate the volume of a cylinder with radius r and length l
	 * 
	 * @param r Radius of the cylinder
	 * @param l Length of the cylinder
	 * @return volume of cylinder
	 */
	public static final double volumeOfACylinder(double r, double l) {
		return 3.1416d*r*r*l;
	}

	/**
	 * \brief Calculate the area of circle with radius r
	 * 
	 * Calculate the area of circle with radius r
	 * 
	 * @param r Radius of the circle
	 * @return area of circle
	 */
	public static final double areaOfACircle(double r) {
		return 3.1416d*r*r;
	}

	/**
	 * \brief Calculate the radius of a sphere with volume v
	 * 
	 * @param v Volume of the sphere
	 * @return Radius of the sphere
	 */
	public static final double radiusOfASphere(double v) {
		return (double) Math.pow(0.23873d*v, 0.33333d);
	}

	/**
	 * \brief Calculate the radius of a sphere with volume v assuming repetition along 3rd dimension
	 * 
	 * Calculate the radius of a sphere with volume v assuming repetition along 3rd dimension
	 * 
	 * @param v Volume of the sphere
	 * @param lZ	Omitted dimension
	 * @return Radius of the sphere
	 */
	public static final double radiusOfASphere(double v, double lZ) {
		return (double) Math.pow(0.23873d*v*lZ, 0.25d);
	}

	/**
	 * \brief Returns the radius of a cylinder with volume v and length l
	 * 
	 * Returns the radius of a cylinder with volume v and length l
	 * 
	 * @param v	Volume of the cylinder
	 * @param l	Length of the cylinder
	 * @return	Radius of the cylinder
	 */
	public static final double radiusOfACylinder(double v, double l) {
		return (double) Math.sqrt(v/(3.1416d*l));
	}

	/**
	 * \brief Calculate the distance between 2 points for points specified in X,Y,Z
	 * 
	 * Calculate the distance between 2 points for points specified in X,Y,Z
	 * 
	 * @param x1	X coordinate of point 1
	 * @param y1	Y coordinate of point 1
	 * @param z1	Z coordinate of point 1
	 * @param x2	X coordinate of point 2
	 * @param y2	Y coordinate of point 2
	 * @param z2	Z coordinate of point 2
	 * @return Distance between the two points
	 */
	public static final double pointDistance(double x1, double y1, double z1, double x2, double y2,
	        double z2) {
		return (double) Math.sqrt(sq(x1-x2)+sq(y1-y2)+sq(z1-z2));
	}

	/**
	 * \brief Calculate the distance between 2 points for points specified as continuous vectors
	 * 
	 * Calculate the distance between 2 points for points specified as continuous vectors
	 * 
	 * @param p1	Point 1
	 * @param p2	Point 2
	 * @return Distance between the two points
	 */
	public static final double pointDistance(ContinuousVector p1, ContinuousVector p2) {
		return pointDistance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
	}

	/**
	 * \brief Perform gamma correction of value v: v^gamma check that v is in the [0,* 1] range
	 * 
	 * Perform gamma correction of value v: v^gamma check that v is in the [0,* 1] range
	 * 
	 * @param v	Value to check
	 * @param gamma	Gamma Constant
	 * @return V corrected by gamma
	 */
	public static double gammaCorrection(double v, double gamma) {
		if ((v<0)|(v>1)) throw new ModelRuntimeException("invalid v for gamma correction, v = "+v);
		return (double) Math.pow(v, gamma);
	}

	/**
	 * \brief Return the maximum among 2 doubles
	 * 
	 * Return the maximum among 2 doubles
	 * 
	 * @param a	Double 1
	 * @param b Double 2
	 * @return the maximum among a and b
	 */
	public static double max(double a, double b) {
		return (a>b ? a : b);
	}

	/**
	 * \brief Return the minimum among 2 doubles
	 * 
	 * Return the minimum among 2 doubles
	 * 
	 * @param a	Double 1
	 * @param b Double 2
	 * @return the minimum among a and b
	 */
	public static double min(double a, double b) {
		return (a<b ? a : b);
	}

	/**
	 * \brief Return the maximum entry in a double array
	 * 
	 * Return the maximum entry in a double array
	 * 
	 * @param a Vector of doubles
	 * @return max double in a
	 */
	public static double max(double[] a) {
		double out = a[0];
		for (int i = 0; i<a.length; i++) {
			out = max(out, a[i]);
		}
		return out;
	}

	/**
	 * \brief Return the average entry in a double array
	 * 
	 * Return the average entry in a double array
	 * 
	 * @param a Vector of doubles
	 * @return average double in a
	 */
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

	/**
	 * \brief Return the minimum entry in a double array
	 * 
	 * Return the minimum entry in a double array
	 * 
	 * @param a Vector of doubles
	 * @return minimum double in a
	 */
	public static double min(double[] a) {
		double out = a[0];
		for (int i = 0; i<a.length; i++) {
			out = Math.min(out, a[i]);
		}
		return out;
	}

	/**
	 * \brief Return the maximum square among 2 doubles
	 * 
	 * Return the maximum square among 2 doubles
	 * 
	 * @param a	Double 1
	 * @param b	Double 2
	 * @return the maximum square among a and b
	 */
	public static double maxSquare(double a, double b) {
		double a2 = a*a;
		double b2 = b*b;
		return (a2>b2 ? a2 : b2);
	}

	/**
	 * \brief Determine if two doubles are the same sign
	 * 
	 * Determine if two doubles are the same sign
	 * 
	 * @param a	Double 1
	 * @param b	Double 2
	 * @return	Boolean noting whether the two are the same sign
	 */
	static public boolean sameSign(double a, double b) {
		return (a*b>=0);
	}

	/**
	 * \brief Return the sum of a vector
	 * 
	 * Return the sum of a vector
	 * 
	 * @param vector	Vector to sum
	 * @return	The sum of that vector
	 */
	public static double sumVector(double vector[]) {
		double out = 0;
		for (int i = 0; i<vector.length; i++) {
			out += vector[i];
		}
		return out;
	}

	/**
	 * \brief Output a double value as a string, in a particular decimal format
	 * 
	 * Output a double value as a string, in a particular decimal format
	 * 
	 * @param value	Value to be output
	 * @param scFormat	The decimal format to use
	 * @return	A string containing that value in the required decimal format
	 */
	public static String toString(double value, boolean scFormat) {
		if (scFormat) return dfSc.format(value);
		else return dfUs.format(value);
	}
	
	/**
	 * \brief Return a uniform distributed random number between 0 and 1
	 * 
	 * Return a uniform distributed random number between 0 and 1
	 * 
	 * @return uniform distributed random number in [0,1]
	 */
	public static double getUniRand() {
		return random.nextDouble();
	}

	/**
	 * \brief Return an integer random number between two set bounds
	 * 
	 * Return an integer random number between two set bounds
	 * 
	 * @param lBound	Lower bound
	 * @param hBound	Upper bound
	 * @return an int between lBound inclusive and hBound exclusive
	 */
	public static int getUniRandInt(int lBound, int hBound) {
		return random.nextInt(hBound-lBound)+lBound;
	}
	
	/**
	 * \brief Return a double random number between two set bounds
	 * 
	 * Return a double random number between two set bounds
	 * 
	 * @param lBound	Lower bound
	 * @param hBound	Upper bound
	 * @return an int between lBound inclusive and hBound exclusive
	 */
	public static double getUniRand(double lBound, double hBound) {
		return random.nextDouble()*(hBound-lBound)+lBound;
	}

	/**
	 * \brief Return a truncated N(0,1) distributed random number. Normal distributed random numbers are truncated at 2*sigma to prevent extreme values.
	 * 
	 * Return a truncated N(0,1) distributed random number. Normal distributed random numbers are truncated at 2*sigma to prevent extreme values.
	 * 
	 * @return truncated N(0,1) distributed random number. 
	 */
	public static double getNormRand() {
		double phi;
		do {
			phi = random.nextGaussian();
		} while (Math.abs(phi)>2);
		return phi;
	}
	
	/**
	 * \brief Return 2 to the power of a uniformly distributed random number in [0,1]
	 * 
	 * Return 2 to the power of a uniformly distributed random number in [0,1]
	 * 
	 * @return 2 to the power of a uniformly distributed random number in [0,1]
	 */
	public static double getExp2Rand() {
		return exp2(getUniRand());
	}
	
	/**
	 * \brief Randomise a value with a gaussian distribution in a range fixed by the ratio sigma
	 * 
	 * Randomise a value with a gaussian distribution in a range fixed by the ratio sigma. If mu = 1 and sigma = .1, the results 
	 * form a truncated gaussian distribution between 0.8 and 1.2 (2*sigma interval)
	 * 
	 * @param mu mean value
	 * @param sigma standard deviation
	 * @return N(mu,sigma)-distributed random value within [mu-2*sigma,mu+2*sigma]
	 */
	public static double deviateFrom(double mu, double sigma) 
	{
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
	
	/**
	 * \brief Calculate the mean of a double array v
	 * 
	 * Calculate the mean of a double array v
	 * 
	 * @param v	The double array
	 * @return	The mean of the double array
	 */
	public static double mean(double [] v) {
		return sumVector(v)/v.length;
	}

	/**
	 * \brief Calculate the standard deviation of a double array v
	 * 
	 * Calculate the standard deviation of a double array v
	 * 
	 * @param v	The double array
	 * @return	The standard deviation of the double array
	 */
	public static double stddev(double [] v) {
		double mean = mean(v);
		double sum = 0.;
		for (int i = 0; i<v.length; i++) {
			sum += sq(v[i] - mean);
		}

		return Math.sqrt(sum/v.length);
	}

	                                      
}