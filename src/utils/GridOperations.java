
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Collection of static methods to deal with spatial grids that are used in the
 * MicroCoSm model. Includes methods for dealing with the vectorized array data
 * structure and for checking consistency of different cubic grids.
 */

/**
 * 
 * @since August 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 */

package utils;

import exceptions.ModelRuntimeException;

public abstract class GridOperations {

	/**
	 * Transforms a 3D array to a 1D array
	 * 
	 * @param originalArray
	 *            3 dimensional array to transform
	 * @return the vectorized array
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
	 * Transforms a 2D array to a 1D array
	 * 
	 * @param originalArray
	 *            2 dimensional array to transform
	 * @return the vectorized array
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
	 * Reshape a vectorized array to 3 dimensional array
	 * 
	 * @param vectorized
	 *            the vectorized array
	 * @param nI,nJ,nK
	 *            the size of the reshaped array
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
	 * Finding the greatest common divisor of two integer numbers
	 * 
	 * @param a,b
	 *            the two numbers
	 * @return the greatest common divisor (int)
	 * @deprecated
	 */
	public static int gcd(int a, int b) {
		if (b == 0)
			return a;
		return gcd(b, a % b);
	}

	/**
	 * Check the consistency of two grid sizes. Two grids of cubic voxels are
	 * consistent if one could be transformed to the other by dividing and then
	 * multipliying the number of voxels in each dimension by an integer number,
	 * respectively.
	 * 
	 * @param a,b
	 *            the grid sizes of the two grids as an array of integer
	 *            {nI,nJ,nK}
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
	 * Add every entry of array b to the corresponding entry in array a
	 * 
	 * @param a,b
	 */
	public static void addTo(double[] a, double[] b) {
		if (a.length != b.length)
			throw new ModelRuntimeException("Trying to add arrays"
					+ " of different sizes:" + a.length + ", " + b.length);
		for (int i = 0; i < a.length; i++)
			a[i] += b[i];
	}

	/**
	 * Subtract every entry of matrix b from the corresponding entry in matrix
	 * a. (Perform a-b)
	 * 
	 * @param a,b
	 */
	public static void subtractFrom(double[] a, double[] b) {
		if (a.length != b.length)
			throw new ModelRuntimeException("Trying to subtract arrays"
					+ " of different sizes:" + a.length + ", " + b.length);
		for (int i = 0; i < a.length; i++)
			a[i] -= b[i];
	}

	/**
	 * @param a
	 * @return the sum of all elements of a
	 */
	public static double computeSum(double[] a) {
		double sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		return sum;
	}

}
