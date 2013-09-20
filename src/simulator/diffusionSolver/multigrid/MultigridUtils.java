/**
 * \package diffusionSolver.multigrid
 * \brief Package of classes used to aid solver calculation for multi-grid scenarios.
 * 
 * Package of classes used to capture the diffusion solvers that can be defined in the protocol file. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.diffusionSolver.multigrid;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import utils.ExtraMath;

import simulator.SoluteGrid;

public abstract class MultigridUtils {

	// boundary layer threshold
	private static final double BLTHRESH  = MultigridSolute.BLTHRESH;
	public static final String  SEPARATOR = " ";

	/**
	 * retuns log2(n - 1) * 0.5
	 * 
	 * @param n
	 * @return order of multigrid
	 * @throws Exception
	 */
	public static int order(Integer n) throws Exception {
		Double out = ExtraMath.log2(n - 1.0);
		if (out%1==0) return out.intValue();
		else throw (new Exception("invalid grid value ("+n+") must be 2^i + 1"));
	}

	/**
	 * \brief Restricts the data in matrix u to a grid one order coarser. Restriction excludes border points.
	 * 
	 * Restricts the data in matrix u to a grid one order coarser. Restriction excludes border points.
	 * 
	 * @param fineGrid finer grid
	 * @param coarsegrid coarser grid
	 */
	public static void restrict2(SoluteGrid fineGrid, SoluteGrid coarsegrid) {

		double[][][] uc = coarsegrid.grid;
		double[][][] u = fineGrid.grid;

		int lc = uc[0][0].length-2;
		int mc = uc[0].length-2;
		int nc = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (lc==1 ? 1.0f/16.0f : 1.0f/32.0f); // pre-compute

		for (k = 1, kc = 1; kc<=lc; kc++, k += 2)
			for (j = 1, jc = 1; jc<=mc; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nc; ic++, i += 2) {
					// 4-connectivity weight
					uc[ic][jc][kc] = 2*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]);
					// 8-connectivity weight
					uc[ic][jc][kc] += u[i+1][j+1][k]+u[i+1][j-1][k]+u[i-1][j+1][k]+u[i-1][j-1][k];
					// 3rd dimension (4-C)
					uc[ic][jc][kc] += 2*(lc==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
					// 3rd dimension (8-C)					
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i-1][j-1][k-1]+u[i-1][j][k-1]+u[i-1][j+1][k-1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j-1][k-1]+u[i][j+1][k-1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i+1][j-1][k-1]+u[i+1][j][k-1]+u[i+1][j+1][k-1]);
					
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i-1][j-1][k+1]+u[i-1][j][k-1]+u[i-1][j+1][k+1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j-1][k+1]+u[i][j+1][k+1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i+1][j-1][k+1]+u[i+1][j][k-1]+u[i+1][j+1][k+1]);
					
					uc[ic][jc][kc] += 4*u[i][j][k];
					uc[ic][jc][kc] *= nfac;
				}

		coarsegrid.refreshBoundary();
	}

	public static void restrict(SoluteGrid fineGrid, SoluteGrid coarseGrid) {

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int lc = uc[0][0].length-2;
		int mc = uc[0].length-2;
		int nc = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (lc==1 ? 1.0f/8.0f : 1.0f/12.0f); // pre-compute

		for (k = 1, kc = 1; kc<=lc; kc++, k += 2)
			for (j = 1, jc = 1; jc<=mc; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nc; ic++, i += 2) {
					// special case for 2D (when lc = 1)
					uc[ic][jc][kc] = u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k];
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
					uc[ic][jc][kc] *= nfac;
					uc[ic][jc][kc] += 0.5f*u[i][j][k];
					if(Double.isNaN(uc[ic][jc][kc])){
						System.out.print(1);
					}
				}
		coarseGrid.refreshBoundary();
	}

	/**
	 * Restricts the data in matrix u to a grid one order coarser. Restriction
	 * excludes border points for points inside the boundary layer, defined by
	 * data in bl. Restriction excludes border points and points outside the
	 * boundary layer (where bl >= 0.5). Points outside boundary layer are
	 * skipped and, therefore, preserve their original value.
	 * 
	 * @param fineGrid	finer grid
	 * @param coarseGrid	coarser grid
	 * @param bl	boundary layer at coarser grid
	 */
	public static void restrictBoundaryLayer2(SoluteGrid fineGrid, SoluteGrid coarseGrid,
	        double[][][] bl) {

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = uc[0][0].length-2;
		int nJ = uc[0].length-2;
		int nI = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (nK==1 ? 1.0f/16.0f : 1.0f/32.0f); // pre-compute

		for (k = 1, kc = 1; kc<=nK; kc++, k += 2)
			for (j = 1, jc = 1; jc<=nJ; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nI; ic++, i += 2)
					if (bl[ic][jc][kc]>=BLTHRESH) {
						// 4-connectivity weight
						uc[ic][jc][kc] = 2*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]);
						// 8-connectivity weight
						uc[ic][jc][kc] += u[i+1][j+1][k]+u[i+1][j-1][k]+u[i-1][j+1][k]+u[i-1][j-1][k];
						// 3rd dimension (4-C)
						uc[ic][jc][kc] += 2*(nK==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
						// 3rd dimension (8-C)					
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i-1][j-1][k-1]+u[i-1][j][k-1]+u[i-1][j+1][k-1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j-1][k-1]+u[i][j+1][k-1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i+1][j-1][k-1]+u[i+1][j][k-1]+u[i+1][j+1][k-1]);
						
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i-1][j-1][k+1]+u[i-1][j][k-1]+u[i-1][j+1][k+1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j-1][k+1]+u[i][j+1][k+1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i+1][j-1][k+1]+u[i+1][j][k-1]+u[i+1][j+1][k+1]);
						
						uc[ic][jc][kc] += 4*u[i][j][k];
						uc[ic][jc][kc] *= nfac;
					}
		coarseGrid.refreshBoundary();
	}

	public static void restrictBoundaryLayer(SoluteGrid fineGrid, SoluteGrid coarseGrid,
	        double[][][] bl) {

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = uc[0][0].length-2;
		int nJ = uc[0].length-2;
		int nI = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (nK==1 ? 1.0f/8.0f : 1.0f/12.0f); // pre-compute

		for (k = 1, kc = 1; kc<=nK; kc++, k += 2)
			for (j = 1, jc = 1; jc<=nJ; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nI; ic++, i += 2)
					if (bl[ic][jc][kc]>=BLTHRESH) {
						// special case for 2D (when lc = 1)
						uc[ic][jc][kc] = u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k];
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
						uc[ic][jc][kc] *= nfac;
						uc[ic][jc][kc] += 0.5f*u[i][j][k];
					}
		coarseGrid.refreshBoundary();
	}	
	/**
	 * Interpolates the data in matrix uc to a grid one order finer for cubic
	 * matrices. Interpolation excludes border points.
	 * 
	 * @param fineGrid finer grid
	 * @param coarsegrid coarser grid
	 */
	static void interpolate(SoluteGrid fineGrid, SoluteGrid coarsegrid) {
		double[][][] uc = coarsegrid.grid;
		double[][][] u = fineGrid.grid;

		int l = u[0][0].length-2;
		int m = u[0].length-2;
		int n = u.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// copy points
		for (kc = 1, k = 1; k<=l; kc++, k += 2) {
			for (jc = 1, j = 1; j<=m; jc++, j += 2) {
				for (ic = 1, i = 1; i<=n; ic++, i += 2)
					u[i][j][k] = uc[ic][jc][kc];
			}
		}
		// interpolate vertically
		for (k = 1; k<=l; k += 2) {
			for (j = 1; j<=m; j += 2) {
				for (i = 2; i<n; i += 2)
					u[i][j][k] = 0.5f*(u[i+1][j][k]+u[i-1][j][k]);
			}
		}
		// interpolate sideways
		for (k = 1; k<=l; k += 2) {
			for (j = 2; j<m; j += 2) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j+1][k]+u[i][j-1][k]);
			}
		}
		for (k = 2; k<l; k += 2) {
			for (j = 1; j<=m; j++) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j][k+1]+u[i][j][k-1]);
			}
		}

		fineGrid.refreshBoundary();
	}

	static void interpolate2(SoluteGrid fineGrid, SoluteGrid coarsegrid) {
		double[][][] uc = coarsegrid.grid;
		double[][][] u = fineGrid.grid;

		int l = u[0][0].length-2;
		int m = u[0].length-2;
		int n = u.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// copy points
		for (kc = 1, k = 1; k<=l; kc++, k += 2) {
			for (jc = 1, j = 1; j<=m; jc++, j += 2) {
				for (ic = 1, i = 1; i<=n; ic++, i += 2)
					u[i][j][k] = uc[ic][jc][kc];
			}
		}
		// interpolate vertically
		for (k = 1; k<=l; k += 2) {
			for (j = 1; j<=m; j += 2) {
				for (i = 2; i<n; i += 2)
					u[i][j][k] = 0.5f*(u[i+1][j][k]+u[i-1][j][k]);
			}
		}
		// interpolate sideways
		for (k = 1; k<=l; k += 2) {
			for (j = 2; j<m; j += 2) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j+1][k]+u[i][j-1][k]);
			}
		}
		for (k = 2; k<l; k += 2) {
			for (j = 1; j<=m; j++) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j][k+1]+u[i][j][k-1]);
			}
		}

		fineGrid.refreshBoundary();
	}
	
	/**
	 * Interpolates the data in matrix uc to a grid one order finer for cubic
	 * matrices for points inside the boundary layer, defined by data in bl.
	 * Interpolation excludes border points and points outside the boundary
	 * layer (where bl >= 0.5). Points outside boundary layer are skipped and,
	 * therefore, preserve their original value.
	 * 
	 * @param fineGrid finer grid
	 * @param coarseGrid coarser grid
	 * @param bl boundary layer at finer grid
	 */
	static void interpolateBoundaryLayer(SoluteGrid fineGrid, SoluteGrid coarseGrid, double[][][] bl) {

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = u[0][0].length-2;
		int nJ = u[0].length-2;
		int nI = u.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// copy points
		for (kc = 1, k = 1; k<=nK; kc++, k += 2) {
			for (jc = 1, j = 1; j<=nJ; jc++, j += 2) {
				for (ic = 1, i = 1; i<=nI; ic++, i += 2)
					if (bl[i][j][k]>=BLTHRESH) u[i][j][k] = uc[ic][jc][kc];
			}
		}
		// interpolate verically
		for (k = 1; k<=nK; k += 2) {
			for (j = 1; j<=nJ; j += 2) {
				for (i = 2; i<nI; i += 2)
					if (bl[i][j][k]>=BLTHRESH) u[i][j][k] = 0.5f*(u[i+1][j][k]+u[i-1][j][k]);
			}
		}
		// interpolate sideways
		for (k = 1; k<=nK; k += 2) {
			for (j = 2; j<nJ; j += 2) {
				for (i = 1; i<=nI; i++)
					if (bl[i][j][k]>=BLTHRESH) u[i][j][k] = 0.5f*(u[i][j+1][k]+u[i][j-1][k]);
			}
		}
		for (k = 2; k<nK; k += 2) {
			for (j = 1; j<=nJ; j++) {
				for (i = 1; i<=nI; i++)
					if (bl[i][j][k]>=BLTHRESH) u[i][j][k] = 0.5f*(u[i][j][k+1]+u[i][j][k-1]);
			}
		}

		fineGrid.refreshBoundary();
	}

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
	 * Set all entries of a boolean matrix to value val
	 * 
	 * @param u
	 * @param val
	 */
	public static void setValues(boolean u[][][], boolean val) {
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
	static void addTo(double a[][][], double b[][][]) {
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] += b[i][j][k];
	}

	/**
	 * Subtract every entry of matrix b to the corresponding entry in matrix a.
	 * 
	 * @param a
	 * @param b
	 */
	static void subtractTo(double a[][][], double b[][][]) {
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
	public static float[][][] subtract(float a[][][], float b[][][]) {
		int l = a.length;
		int m = a[0].length;
		int n = a[0][0].length;
		float[][][] c = new float[l][m][n];
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
	public static float min(float a[][][]) {
		float min = a[0][0][0];
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
	public static float max(float a[][][]) {
		float max = a[0][0][0];
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
	public static float computeNorm(double[][][] a) {
		float norm = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					norm += ExtraMath.sq(a[i][j][k]);
		return (float) Math.sqrt(norm);
	}

	/**
	 * @param a
	 * @return the sum of all elements of a
	 */
	public static float computeSum(float[][][] a) {
		float sum = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					sum += a[i][j][k];
		return sum;
	}

	/**
	 * Return values in a matrix (excluding boundaries) as a formatted string
	 * 
	 * @param matrix to output as string
	 * @return string output
	 */
	public static String coreMatrixToString(float[][][] matrix) {
		int n = matrix.length-2;
		int m = matrix[0].length-2;
		int l = matrix[0][0].length-2;
		StringBuffer out = new StringBuffer();
		for (int k = 1; k<=l; k++) {
			for (int i = n; i>=1; i--) {
				for (int j = 1; j<=m; j++) {
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
	 * Return values in a matrix (excluding boundaries) as a formatted string.
	 * This method is used for boolean matrices. Values in output are 1 (for
	 * true) or 0 (for false)
	 * 
	 * @param matrix to output as string
	 * @return string output
	 */
	public static String coreMatrixToString(boolean[][][] matrix) {
		int n = matrix.length-2;
		int m = matrix[0].length-2;
		int l = matrix[0][0].length-2;
		StringBuffer out = new StringBuffer();
		for (int k = 1; k<=l; k++) {
			for (int i = n; i>=1; i--) {
				for (int j = 1; j<=m; j++) {
					out.append(matrix[i][j][k] ? 1 : 0);
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
	 * Write the full matrix to a string
	 * 
	 * @param matrix
	 * @return a string with the matrix (space separated values)
	 */
	public static String matrixToString(float[][][] matrix) {
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
	 * Create a 2D graphics
	 * 
	 * @param fileName the file to parse
	 * @return 2D matrix
	 */
	public static float[][] readSquareMatrixFromFile(String fileName) {
		String line;
		String[] tokens;
		float[][] dataRead = { { 0 } }; // initialize
		int lineCount = 0;
		int i = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while ((line = br.readLine())!=null) {
				tokens = line.split(SEPARATOR);
				int l = tokens.length;
				// in first iteration initialize the 2D matrix
				if (dataRead.length==1) {
					// create the matrix and
					lineCount = l;
					dataRead = new float[lineCount][lineCount];
					i = lineCount-1;
				} else if ((lineCount!=l)&(i>=0)) {

				} else if (i<-1) {

				}
				if (i>=0) {
					// parse the data in the line into a matrix
					for (int j = 0; j<lineCount; j++) {
						dataRead[i][j] = Float.parseFloat(tokens[j]);
					}
				}
				// decrement the line number counter
				i--;
			}
		} catch (FileNotFoundException e) {

		} catch (IOException e) {

		}
		return dataRead;
	}

}