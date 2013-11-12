/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of
 * running an iDynoMiCS Simulation.
 * 
 * Package of classes that perform utility functions in the process of running
 * an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by
 * the CeCILL license under French law and abides by the rules of distribution
 * of free software.  You can use, modify and/ or redistribute iDynoMiCS under
 * the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the
 * following URL "http://www.cecill.info".
 */
package utils;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Random;

/**
 * \brief Abstract class with some extra useful math functions.
 * 
 * Contents:
 * 		Simple calculations
 * 		Shapes
 * 		Arrays of doubles
 * 		Dealing with signs
 * 		Dealing with strings
 * 		Random number generation
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA)
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
 * @author Robert Clegg (rjc096@bham.ac.uk), University of Birmingham, UK
 */
public final class ExtraMath
{
	/**
	 * \brief One option for writing decimals to screen.
	 * 
	 * This always has 3 digits before the decimal point, and adjusts the
	 * scientific exponent accordingly.
	 */
	public static java.text.DecimalFormat dfSc = new DecimalFormat("000.###E0");
	
	/**
	 * \brief Second option for writing decimals to screen.
	 * 
	 * This always has 2 digits after the decimal point, and will round any
	 * smaller decimals.
	 */
	public static java.text.DecimalFormat dfUs = new DecimalFormat("########.##");
	
	/**
	 * Random number generator
	 */
	public static Random random;
	
	/* ----------------------- Simple calculations ------------------------ */
	
	/**
	 * \brief Computes the logarithm of base 2.
	 * 
	 * If x is non-positive Double.NaN will be returned.
	 * 
	 * @param x The Integer to take the logarithm of.
	 * @return Double the logarithm base 2 of x
	 */
	public static final Double log2(Double x)
	{
		return (Math.log(x)/Math.log(2));
	}
	
	/**
	 * \brief Square an Integer number.
	 * 
	 * @param x The Integer to square.
	 * @return The Integer square of x.
	 */
	public static final Integer sq(Integer x)
	{
		return x*x;
	}
	
	/**
	 * \brief Square a Double number.
	 * 
	 * @param x The Double to square.
	 * @return The Double square of x.
	 */
	public static final Double sq(Double x) 
	{
		return x*x;
	}
	
	/**
	 * \brief Cube an Integer number.
	 * 
	 * @param x The Integer to cube.
	 * @return The Integer cube of x.
	 */
	public static final Integer cube(Integer x)
	{
		return x*x*x;
	}
	
	/**
	 * \brief Cube a Double number.
	 * 
	 * @param x The Double to cube.
	 * @return The Double cube of x.
	 */
	public static final Double cube(Double x)
	{
		return x*x*x;
	}
	
	/**
	 * \brief Find the real cube root of a double number.
	 * 
	 * @param x The Double to take the cube root of.
	 * @return The Double real cube root of x.
	 */
	public static final Double cubeRoot(Double x)
	{
		return Math.pow(x, (1.0/3.0));
	}
	
	/**
	 * \brief Calculate 2 to the power of x where x is an integer.
	 * 
	 * Returns 1 if x is less than zero.
	 * 
	 * @param x The exponent
	 * @return 2^x
	 */
	public static final Integer exp2(Integer x)
	{
		Integer out = 1;
		for (int i = 0; i < x; i++)
			out *= 2;
		return out;
	}
	
	/**
	 * \brief Calculate 2 to the power of x where x is a long integer.
	 * 
	 * Returns 1 if x is less than zero.
	 * 
	 * @param x
	 * @return 2^x
	 */
	public static final long exp2long(int x)
	{
		long out = 1;
		for (int i = 0; i< x; i++)
			out *= 2;
		return out;
	}
	
	/**
	 * \brief Calculate 2 to the power of x where x is a double.
	 * 
	 * @param x The exponent
	 * @return 2^x
	 */
	public static final Double exp2(Double x)
	{
		return Math.pow(2, x);
	}
	
	/*  ----------------------------- Shapes  ----------------------------- */
	
	/**
	 * \brief Calculate the area of circle with radius r.
	 * 
	 * The formula for this is pi r^2.
	 * 
	 * @param radius Radius of the circle
	 * @return area of circle
	 */
	public static final Double areaOfACircle(Double radius)
	{
		return Math.PI * sq(radius);
	}
	
	/**
	 * \brief Calculate the volume of a cylinder with radius r and length l.
	 * 
	 * The formula for this is pi * r^2 * l.
	 * 
	 * @param radius Radius of the cylinder
	 * @param length Length of the cylinder
	 * @return volume of cylinder
	 */
	public static final Double volumeOfACylinder(Double radius, Double length)
	{
		return areaOfACircle(radius) * length;
	}
	
	/**
	 * \brief Calculate the radius of a circle with area a.
	 * 
	 * The formula for this is ( a/pi )^(1/2)
	 * 
	 * @param area Area of the circle
	 * @return Radius of the circle
	 */
	public static final Double radiusOfACircle(Double area)
	{
		return Math.sqrt(area / Math.PI);
	}
	
	/**
	 * \brief Calculate the radius of a cylinder with volume v and length l.
	 * 
	 * This is calculated from the area of the cross-section: v/l
	 * 
	 * @param volume Volume of the cylinder
	 * @param length Length of the cylinder
	 * @return Radius of the cylinder
	 */
	public static final Double radiusOfACylinder(Double volume, Double length)
	{
		return radiusOfACircle(volume/length);
	}
	
	/**
	 * \brief Calculate the volume of a sphere with radius r.
	 * 
	 * The formula for this is 4/3 * pi * r^3.
	 * 
	 * @param radius Radius of the sphere
	 * @return volume of sphere
	 */
	public static final Double volumeOfASphere(Double radius)
	{
		return (4.0/3.0) * Math.PI * cube(radius);
	}
	
	/**
	 * \brief Calculate the radius of a sphere with volume v.
	 * 
	 * The formula for this is ( (v*3)/(4*pi) )^(1/3)
	 * 
	 * @param volume Volume of the sphere
	 * @return Radius of the sphere
	 */
	public static final Double radiusOfASphere(Double volume)
	{
		return cubeRoot(volume*0.75/Math.PI);
	}
	
	/*  ------------------------ Arrays of doubles  ------------------------ */
	
	/**
	 * \brief Set up a new Double[] that is initialised with zero's.
	 * 
	 * Note that initialising a Double[] will normally fill it with null's.
	 * 
	 * @param length Integer length of the array requested.
	 * @return Double[] array of zero's
	 */
	public static Double[] newDoubleArray(Integer length)
	{
		return Collections.nCopies(length, 0.0).toArray(new Double[0]);
	}
	
	/**
	 * \brief Return the maximum entry in a Double array.
	 * 
	 * @param array Array of Doubles.
	 * @return Greatest Double in a.
	 */
	public static Double max(Double[] array)
	{
		Double out = array[0];
		for (Double value : array)
			out = Math.max(out, value);
		return out;
	}
	
	/**
	 * \brief Return the minimum entry in a Double array.
	 * 
	 * @param array Array of Doubles.
	 * @return Least Double in a.
	 */
	public static Double min(Double[] array)
	{
		Double out = array[0];
		for (Double value : array)
			out = Math.min(out, value);
		return out;
	}
	
	/**
	 * \brief Determine the sum of an array of doubles.
	 * 
	 * Takes care to exclude any infinites or NaN's.
	 * 
	 * @param array Array of doubles
	 * @return
	 */
	public static Double sum(Double[] array)
	{
		Double sum = 0.0;
		for (Double value : array)
			if(!Double.isInfinite(value) && !Double.isNaN(value))
				sum += value;
		return sum;
	}
	
	/**
	 * \brief Return the mean average entry in a Double array.
	 * 
	 * Removes infinite/NaN values from the sum and the number of entries.
	 * 
	 * @param array Vector of Doubles.
	 * @return Average Double in array.
	 */
	public static Double mean(Double[] array)
	{
		Double out = 0.0;
		Double n = 0.0;
		for (Double value : array)
			if(!Double.isInfinite(value) && !Double.isNaN(value))
			{
					out += value;
					n++;
			}
		if (n == 0.0)
		{
			LogFile.writeLogAlways("WARNING! ExtraMath.mean(): array of length "+
												array.length+" has no valid entries");
			return 0.0;
		}
		return (out/n);
	}
	
	/**
	 * \brief Return the mean standard deviation of a Double array.
	 * 
	 * Removes infinite/NaN values from the sum and the number of entries.
	 * 
	 * @param array Vector of Doubles.
	 * @param fromSample Boolean denoting whether to divide by n-1 (True) or n (False) 
	 * @return Standard deviation of a Double array
	 */
	public static Double stddev(Double[] array, Boolean fromSample)
	{
		Double mean = mean(array);
		Double sum = 0.0;
		Double n = 0.0;
		
		for (Double value : array)
			if(!Double.isInfinite(value) && !Double.isNaN(value))
			{
				sum += sq(value - mean);
				n++;
			}
		
		// check the array contains valid entries before trying to divide by zero
		if (n == 0.0)
		{
			LogFile.writeLogAlways("WARNING! ExtraMath.stddev(): array of length "+
												array.length+" has no valid entries");
			return 0.0;
		}
		
		// If this is from a sample we divide by (n-1), not n
		if ((fromSample) && (n>1))
			n--;
		return Math.sqrt(sum/n);
	}
	
	/*  ----------------------- Dealing with signs  ----------------------- */
	
	/**
	 * \brief Unequivocally determine the sign of a double number. 
	 * 
	 * Copied from 
	 * http://stackoverflow.com/questions/3994531/how-to-determine-if-a-number-is-positive-or-negative-in-java
	 * on 7 August 2013
	 * 
	 * @param f Double to be inspected
	 * @return Integer value with the sign of f: -1, 0, or +1
	 */
	public static Integer sign(Double f)
	{
	    if (f == 0.0) return 0;
	    f *= Double.POSITIVE_INFINITY;
	    if (f == Double.POSITIVE_INFINITY) return +1;
	    if (f == Double.NEGATIVE_INFINITY) return -1;

	    throw new IllegalArgumentException("Unfathomed double");
	}
	
	/**
	 * \brief Determine if two Doubles are the same sign.
	 * 
	 * Note that this is true if if either (or both) of the arguments is zero.
	 * 
	 * @param a	Double 1
	 * @param b	Double 2
	 * @return	Boolean noting whether the two are the same sign.
	 */
	static public Boolean sameSign(Double a, Double b)
	{
		return (sign(a)*sign(b) >= 0);
	}
	
	/**
	 * \brief Output a Double value as a string, in a particular decimal format.
	 * 
	 * If true, use dfSc; if false, use dfUs. 
	 * 
	 * @param value	Double value to be output.
	 * @param scFormat	The decimal format to use.
	 * @return	A string containing that value in the required decimal format.
	 */
	public static String toString(Double value, Boolean scFormat)
	{
		return (scFormat) ? dfSc.format(value) : dfUs.format(value); 
	}
	
	/**
	 * \brief Searches for a substring within a main string, and returns a double immediately after if it exists.
	 * 
	 * Note that 1.0 will be returned if the substring is not found, the substring is at the very end 
	 * of the main string, or there is no double immediately after.
	 * 
	 * @param mainString The string within which the search will be made
	 * @param subString The substring being searched for
	 * @return The double immediately after subString, if found. If not found, 1.0
	 */
	public static Double doubleAfterSubstring(String mainString, String subString)
	{
		Double out = 1.0;
		if (mainString.contains(subString))
		{
			int startIndex = mainString.indexOf(subString) + subString.length();
			int endIndex = startIndex + 1;
			int maxIndex = mainString.length();
			String potential;
			
			while ((endIndex < maxIndex) && 
					(isNumeric(mainString.substring(startIndex, endIndex+1))))
				endIndex++;
			
			potential = mainString.substring(startIndex, endIndex);
			
			if (isNumeric(potential))
				out = Double.parseDouble(potential); 
		}
		return out;
	}
	

	/**
	 * \brief Checks if the supplied String can be safely parsed as a Double.
	 * 
	 * @param str The string to be tested.
	 * @return True or False depending on the outcome of the test.
	 */
	public static Boolean isNumeric(String str)  
	{  
	  try  
	  {  
	    @SuppressWarnings("unused")
		double d = Double.parseDouble(str);  
	  }  
	  catch(NumberFormatException nfe)
	  {  
	    return false;
	  }  
	  return true;  
	}
	
	/*  -------------------- Random number generation  -------------------- */
	
	/**
	 * \brief Return a uniformly distributed random number between 0 and 1.
	 * 
	 * Lower bound (0) is inclusive, upper bound (1) is exclusive.
	 * 
	 * @return A uniformly distributed random number in [0,1).
	 */
	public static Double getUniRandDbl() {
		return random.nextDouble();
	}
	
	/**
	 * \brief Return a double random number between two set bounds.
	 * 
	 * @param lBound	Lower bound (inclusive).
	 * @param uBound	Upper bound (exclusive).
	 * @return A uniformly distributed random double number in [lBound, uBound).
	 */
	public static Double getUniRandDbl(Double lBound, Double uBound)
	{
		return getUniRandDbl()*(uBound-lBound)+lBound;
	}
	
	/**
	 * \brief Return 2 to the power of a uniformly distributed random number
	 * in [0,1).
	 * 
	 * @return 2 to the power of a uniformly distributed random number in [0,1).
	 */
	public static Double getExp2Rand()
	{
		return exp2(getUniRandDbl());
	}
	
	/**
	 * \brief Return an integer random number less than the upper bound supplied.
	 * 
	 * @param uBound	Upper bound (exclusive).
	 * @return A uniformly distributed random integer number in [0, uBound).
	 */
	public static Integer getUniRandInt(Integer uBound)
	{
		return random.nextInt(uBound);
	}
	
	/**
	 * \brief Return an integer random number between two set bounds.
	 * 
	 * @param lBound	Lower bound (inclusive).
	 * @param uBound	Upper bound (exclusive).
	 * @return A uniformly distributed random integer number in [lBound, uBound).
	 */
	public static Integer getUniRandInt(Integer lBound, Integer hBound)
	{
		return getUniRandInt(hBound-lBound)+lBound;
	}
	
	/**
	 * @param lBound
	 * @param hBound
	 * @return a double between lBound inclusive and hBound exclusive
	 */
	public static Double getUniRand(double lBound, double hBound) {
		return random.nextDouble()*(hBound-lBound)+lBound;
	}

	/**
	 * \brief Return a truncated N(0,1) distributed random number.
	 * 
	 * Normal distributed random numbers are truncated at 2*sigma to prevent
	 * extreme values.
	 * 
	 * @return Truncated N(0,1) distributed random number. 
	 */
	public static Double getNormRand()
	{
		Double phi;
		do {
			phi = random.nextGaussian();
		} while (Math.abs(phi)>2);
		return phi;
	}
	
	/**
	 * \brief Randomise a value with a normal distribution in a range fixed by
	 * the Coefficient of Variation (CV).
	 * 
	 * Randomise a value mu with a normal (Gaussian) distribution in a range
	 * fixed by CV. This is different from deviateFromSD()!
	 * The result will be the same sign (+/-) as mu.
	 * 
	 * E.g. If mu = 1 and cv = .1, the results form a truncated normal
	 * 			distribution between 0.8 and 1.2
	 * 		If mu = -3 and cv = .05 the results form a truncated normal
	 * 			distribution between -3.3 and -2.7
	 * 
	 * @param mu Mean value.
	 * @param cv Coefficient of Variation.
	 * @return N(mu, cv)-distributed random value within
	 * [mu*(1-2*cv), mu*(1+2*cv)]
	 */
	public static Double deviateFromCV(Double mu, Double cv) 
	{
		// No point going further if either is zero 
		if (mu == 0.0 || cv == 0.0) return mu;
		
		Double result, deviation;
		do {
			deviation = getNormRand();
			// Compute the new value
			result = mu*(1.0+cv*deviation);
		} while (!sameSign(result, mu));
		
		return result;
	}
	
	/**
	 * \brief Randomise a value with a normal distribution in a range fixed by 
	 * the Standard Deviation (SD).
	 * 
	 * Randomise a value mu with a normal (Gaussian) distribution in a range
	 * fixed by SD. This is different from deviateFromCV()!
	 * The result will be the same sign (+/-) as mu.
	 * 
	 * E.g. If mu = 1 and sigma = .1, the results form a truncated normal
	 * 			distribution between 0.8 and 1.2
	 * 		If mu = -3 and sigma = .05 the results form a truncated normal
	 * 			distribution between -3.1 and -2.9
	 * 
	 * @param mu Mean value
	 * @param sigma Standard Deviation
	 * @return N(mu,sigma)-distributed random value within
	 * [mu-2*sigma, mu+2*sigma]
	 */
	public static Double deviateFromSD(Double mu, Double sigma) 
	{
		// No point going further if the standard deviation is zero 
		if (sigma == 0.0) return mu;
		
		Double result, deviation;
		do {
			deviation = getNormRand();
			// Compute the new value
			result = mu + (sigma*deviation);
		} while (!sameSign(result, mu));
		
		return result;
	}                           
}
