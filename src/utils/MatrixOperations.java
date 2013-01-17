/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *______________________________________________________
 * Implements static utility functions for used in multigrid method.
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */

package utils;

public abstract class MatrixOperations {

	static final String SEPARATOR = "\t";

	/**
     * Set all entries of a matrix to value val
     * 
     * @param u
     * @param val
     */
	public static void setValues(double u[][][], double val) {
		for (int i = 0; i<u.length; i++)
			for (int j = 0; j<u[i].length; j++)
				for (int k = 0; k<u[i][j].length; k++)
					u[i][j][k] = val;
	}

	/**
     * Add every entry of matrix b to the corresponding entry in matrix a
     * 
     * @param a
     * @param b
     */
	public static void addTo(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] += b[i][j][k];
	}

	/**
     * Add every entry of matrix b to the corresponding entry in matrix a
     * 
     * @param a
     * @param b
     */
	public static void muliplyBy(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] *= b[i][j][k];
	}

	public static void muliplyBy(double a[][][], double b) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] *= b;
	}

	/**
     * Subtract every entry of matrix b from the corresponding entry in matrix
     * a.
     * 
     * @param a
     * @param b
     */
	public static void subtractFrom(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] -= b[i][j][k];
	}

	/**
     * Create matrix c = a - b
     * 
     * @param a
     * @param b
     * @return c = a-b
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
     * Find minimum value in a 3D matrix
     * 
     * @param a
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
     * Find maximum value in a 3D matrix
     * 
     * @param a
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
     * compute the norm of matrix (exceptuating padding)
     * 
     * @param a
     * @return the norm of the matrix
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
		//return (double) Math.sqrt(norm/(a.length-1)/(a[1].length-1)/(a[1][1].length-1));
		return (double) Math.sqrt(norm);
	}
	
	public static double computeAverage(double[][][] a) {
		return computeSum(a)/(a.length-2)/(a[0].length-2)/(a[0][0].length-2);
	}

	/**
     * @param a
     * @return the sum of all elements of a grid
     * padding excluded
     */
	public static double computeSum(double[][][] a) {
		double sum = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					sum += a[i][j][k];
		return sum;
	}
	
	/*//sonia:chemostat 19.02.2010
	public static double computeSumChemo(double[][][] a) {
		double sum = 0;
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					sum += a[i][j][k];
		return sum;
	}*/

	public static double computeSumP(double[][][] a) {
		double sum = 0;
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					sum += a[i][j][k];
		return sum;
	}

	/**
     * Write the full matrix to a string
     * 
     * @param matrix
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
     * Copy values from an array src to array dest
     * 
     * @param dest write values to this array
     * @param src copy values from this array
     */
	public static void copyValuesTo(double dest[][][], double src[][][]) {
		for (int i = 0; i<dest.length; i++)
			for (int j = 0; j<dest[i].length; j++)
				for (int k = 0; k<dest[i][j].length; k++)
					dest[i][j][k] = src[i][j][k];
	}

}