/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import exceptions.ModelRuntimeException;

/**
 * \brief Collection of static methods to deal with spatial grids that are used in the MicroCoSm model. 
 * 
 * Collection of static methods to deal with spatial grids that are used in the MicroCoSm model. Includes methods for dealing with 
 * the vectorized array data structure and for checking consistency of different cubic grids. 
 * 
 * @author Andreas Dotsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 *
 */
public abstract class GridOperations {

	/**
	 * \brief Transforms a 3D array to a 1D array
	 * 
	 * Transforms a 3D array to a 1D array
	 * 
	 * @param originalArray	3 dimensional array to transform
	 * @return Vectorized array in one dimension
	 * @deprecated
	 */
	public static double[] vectorizeArray(double[][][] originalArray) {
		double[] transformed = new double[originalArray.length
				* originalArray[0].length * originalArray[0][0].length];
		for (int i = 0; i < originalArray.length; i++) {
			for (int j = 0; j < originalArray[0].length; j++) {
				for (int k = 0; k < originalArray[0][0].length; k++) {
					transformed[i + originalArray.length * j
							+ originalArray.length * originalArray[0].length
							* k] = originalArray[i][j][k];
				}
			}
		}
		return transformed;
	}

	/**
	 * \brief Transforms a 2D array to a 1D array
	 * 
	 * Transforms a 2D array to a 1D array
	 * 
	 * @param originalArray	2 dimensional array to transform
	 * @return Vectorized array in one dimension
	 * @deprecated
	 */
	public static double[] vectorizeArray(double[][] originalArray) {
		double[] transformed = new double[originalArray.length
				* originalArray[0].length];
		for (int i = 0; i < originalArray.length; i++) {
			for (int j = 0; j < originalArray[0].length; j++) {
				transformed[i + originalArray.length * j] = originalArray[i][j];
			}
		}
		return transformed;
	}

	/**
	 * \brief Reshape a vectorized array to 3 dimensional array
	 * 
	 * Reshape a vectorized array to 3 dimensional array
	 * 
	 * @param vectorized	the vectorized array
	 * @param nI	Array size in I direction
	 * @param nJ	Array size in J direction
	 * @param nK	Array size in K direction
	 * @return the reshaped array
	 * @deprecated
	 */
	public static double[][][] reshapeArray(double[] vectorized, int nI,
			int nJ, int nK) {
		if (nI * nJ * nK != vectorized.length)
			throw new ModelRuntimeException(
					"Trying to reshape an array of length" + vectorized.length
							+ " to size " + nI + "*" + nJ + "*" + nK);
		double[][][] reshaped = new double[nI][nJ][nK];
		int ii = 0;
		for (int i = 0; i < nI; i++) {
			for (int j = 0; j < nJ; j++) {
				for (int k = 0; k < nK; k++) {
					reshaped[i][j][k] = vectorized[ii];
					ii++;
				}
			}
		}
		return reshaped;
	}

	/**
	 * \brief Find the greatest common divisor of two integer numbers
	 * 
	 * Find the greatest common divisor of two integer numbers
	 * 
	 * @param a	Integer one
	 * @param b	Integer two
	 * @return the greatest common divisor (int)
	 * @deprecated
	 */
	public static int gcd(int a, int b) {
		if (b == 0)
			return a;
		return gcd(b, a % b);
	}

	/**
	 * \brief Check the consistency of two grid sizes.
	 * 
	 * Check the consistency of two grid sizes. Two grids of cubic voxels are consistent if one could be transformed to the other 
	 * by dividing and then multipliying the number of voxels in each dimension by an integer number, respectively.
	 * 
	 * @param a	Grid sizes of first grid as an array of integers {nI,nJ,nK}
	 * @param b	Grid sizes of second grid as an array of integers {nI,nJ,nK}
	 * @return true if grids are consistent
	 */
	public static boolean isConsistent(int[] a, int[] b) {
		if (a.length != b.length)
			throw new ModelRuntimeException("Non-matching grid dimensions!");

		// check size ratios in all dimensions
		float[] ab = new float[a.length];
		ab[0] = (float) a[0] / (float) b[0];
		for (int i = 1; i < a.length; i++) {
			ab[i] = (float) a[i] / (float) b[i];
			if (ab[i] != ab[i - 1])
				return false;
		}

		// find the greatest common divisor of the grid sizes of grid a
		int d = a[0];
		for (int i = 1; i < a.length; i++) {
			d = gcd(d, a[i]);
		}

		// check, if b can be transformed from a coarsened grid a
		for (int i = 0; i < b.length; i++) {
			if (((b[i] % (a[i] / d)) != 0)) {
				return false;
			}
		}

		return true;
	}

	/**
	 * \brief Add every entry of array b to the corresponding entry in array a
	 * 
	 * Add every entry of array b to the corresponding entry in array a
	 * 
	 * @param a	Target array
	 * @param b	Array of values being added to a
	 */
	public static void addTo(double[] a, double[] b) {
		if (a.length != b.length)
			throw new ModelRuntimeException("Trying to add arrays"
					+ " of different sizes:" + a.length + ", " + b.length);
		for (int i = 0; i < a.length; i++)
			a[i] += b[i];
	}

	/**
	 * \brief Subtract every entry of matrix b from the corresponding entry in matrix a. (Perform a-b)
	 * 
	 * Subtract every entry of matrix b from the corresponding entry in matrix a. (Perform a-b)
	 * 
	 * @param a	Target array
	 * @param b	Array of values being subtracted from a
	 */
	public static void subtractFrom(double[] a, double[] b) {
		if (a.length != b.length)
			throw new ModelRuntimeException("Trying to subtract arrays"
					+ " of different sizes:" + a.length + ", " + b.length);
		for (int i = 0; i < a.length; i++)
			a[i] -= b[i];
	}

	/**
	 * \brief Return the sum of all elements in a double array
	 * 
	 * Return the sum of all elements in a double array
	 * 
	 * @param a	Double array to sum
	 * @return the sum of all elements of a
	 */
	public static double computeSum(double[] a) {
		double sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		return sum;
	}

}
