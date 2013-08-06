/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

/**
 * \brief Implements static utility functions for used in multigrid method.
 * 
 * Implements static utility functions for used in multigrid method.
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public abstract class MatrixOperations 
{
	/**
	 * Used where this matrix is converted to a string, to separate the elements
	 */
	static final String SEPARATOR = "\t";

	/**
     * \brief Set all entries of a specified matrix to a provided value
     * 
     * Set all entries of a specified matrix to a provided value
     * 
     * @param u	The matrix for which all values should be set
     * @param val	The value for which all elements should be set
     */
	public static void setValues(double u[][][], double val) {
		for (int i = 0; i<u.length; i++)
			for (int j = 0; j<u[i].length; j++)
				for (int k = 0; k<u[i][j].length; k++)
					u[i][j][k] = val;
	}

	/**
     * \brief Add every entry of matrix b to the corresponding entry in matrix a
     * 
     * Add every entry of matrix b to the corresponding entry in matrix a
     * 
     * @param a	Matrix whos value is being increased
     * @param b	Matrix of values to add to the matrix above
     */
	public static void addTo(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] += b[i][j][k];
	}

	/**
     * \brief Multiply every entry of matrix a by the corresponding entry in matrix b
     * 
     * Multiply every entry of matrix a to the corresponding entry in matrix b
     * 
     * @param a	Matrix whos value is being increased
     * @param b	Matrix of values to multiply to the matrix above
     */
	public static void muliplyBy(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] *= b[i][j][k];
	}

	/**
     * \brief Multiply every entry of matrix a by a specified value b
     * 
     * Multiply every entry of matrix a by a specified value b
     * 
     * @param a	Matrix whos value is being increased
     * @param b	Double value by which every entry in the matrix should be multiplied
     */
	public static void muliplyBy(double a[][][], double b) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] *= b;
	}

	/**
     * \brief Subtract every entry of matrix b from the corresponding entry in matrix a.
     * 
     * Subtract every entry of matrix b from the corresponding entry in matrix a.
     * 
     * @param a	Matrix whos value is being decreased
     * @param b	Matrix containing the values by which to decrease the matrix above
     */
	public static void subtractFrom(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] -= b[i][j][k];
	}

	/**
     * \brief Creates a matrix c by subtracting matrix b from matrix a
     * 
     * @param a	Matrix whos value is being decreased
     * @param b	Matrix containing the values by which to decrease the matrix above
     * @return c = a-b	The matrix created by this subtraction
     */
	public static double[][][] subtract(double a[][][], double b[][][]) {
		int l = a.length;
		int m = a[0].length;
		int n = a[0][0].length;
		double[][][] c = new double[l][m][n];
		for (int i = 0; i<l; i++)
			for (int j = 0; j<m; j++)
				for (int k = 0; k<n; k++)
					c[i][j][k] = a[i][j][k]-b[i][j][k];
		return c;
	}

	/**
     * \brief Find minimum value in a 3D matrix
     * 
     * Find minimum value in a 3D matrix
     * 
     * @param a	The matrix to find the minimum value within
     * @return the minimum value in the matrix
     */
	public static double min(double a[][][]) {
		double min = a[0][0][0];
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					min = (a[i][j][k]<min ? a[i][j][k] : min);
		return min;
	}

	/**
     * \brief Find maximum value in a 3D matrix
     * 
     * Find maximum value in a 3D matrix
     * 
     * @param a	The matrix to find the maximum value within
     * @return the maximum value in the matrix
     */
	public static double max(double a[][][]) {
		double max = a[0][0][0];
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					max = (a[i][j][k]>max ? a[i][j][k] : max);
		return max;
	}

	/**
     * \brief Compute the norm of a given matrix (exceptuating padding)
     * 
     * Compute the norm of a given matrix (exceptuating padding)
     * 
     * @param a	The matrix for which the norm should be computed
     * @return Double value that is the norm of the matrix
     */
	public static double computeNorm(double[][][] a) {
		double norm = 0;
		for (int i = 1; i<a.length-1; i++) {
			for (int j = 1; j<a[i].length-1; j++) {
				for (int k = 1; k<a[i][j].length-1; k++) {
					norm += ExtraMath.sq(a[i][j][k]);
				}
			}
		}
		return (double) Math.sqrt(norm);
	}
	
	/**
     * \brief Compute the average of a given matrix (exceptuating padding)
     * 
     * Compute the average of a given matrix (exceptuating padding)
     * 
     * @param a	The matrix for which the average should be computed
     * @return Double value that is the average of the matrix
     */
	public static double computeAverage(double[][][] a) {
		return computeSum(a)/(a.length-2)/(a[0].length-2)/(a[0][0].length-2);
	}

	/**
     * \brief Compute and return the sum of all elements in the grid, padding excluded
     * 
     * Compute and return the sum of all elements in the grid, padding excluded
     * 
     * @param a	The matrix for which the sum should be computed
     * @return the sum of all elements of a grid
     */
	public static double computeSum(double[][][] a) {
		double sum = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					sum += a[i][j][k];
		return sum;
	}
	
	
	public static double computeSumP2(double[][][] a) {
		double sum = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
				{
					if(a[i][j][k]!=0)
					{
						System.out.println(a[i][j][k]);
					}
				
					sum += a[i][j][k];
				}
		return sum;
	}
	
	/**
     * \brief Compute and return the sum of all elements in the grid, padding included
     *
     * Compute and return the sum of all elements in the grid, padding included
     * 
     * @param a	The matrix for which the sum should be computed
     * @return the sum of all elements of a grid
     */
	public static double computeSumP(double[][][] a) {
		double sum = 0;
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					sum += a[i][j][k];
		return sum;
	}

	/**
     * \brief Write the full matrix to a string
     * 
     * Write the full matrix to a string
     * 
     * @param matrix	The matrix which should be written to a string
     * @return a string with the matrix (space separated values)
     */
	public static String matrixToString(double[][][] matrix) {
		StringBuffer out = new StringBuffer();
		for (int k = 0; k<matrix[0][0].length; k++) {
			for (int i = matrix.length-1; i>=0; i--) {
				for (int j = 0; j<matrix[0].length; j++) {
					out.append(matrix[i][j][k]);
					// change here for format (presently space separated values
					out.append(SEPARATOR);
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	/**
     * \brief Copy values from an array source to an array destination
     * 
     * Copy values from an array source to an array destination
     * 
     * @param dest[][][] Array into which the values should be written
     * @param src[][][] Array into which the values should be written
     */
	public static void copyValuesTo(double dest[][][], double src[][][]) {
		for (int i = 0; i<dest.length; i++)
			for (int j = 0; j<dest[i].length; j++)
				for (int k = 0; k<dest[i][j].length; k++)
					dest[i][j][k] = src[i][j][k];
	}

}