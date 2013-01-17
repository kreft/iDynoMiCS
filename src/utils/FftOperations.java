/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 * ______________________________________________________
 * Tool for performing fast fourier transfrom (FFT) on grids of Complex and real
 * numbers. Algorithms for FFT based on "Numerical Recipes in C, 2nd ed."
 */

/**
 * 
 * @since August 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 */


package utils;


public class FftOperations {

	// TODO: check orginalData transformedData, is it the same because of
    // references?

	public static Complex[] transform(double[] originalData) {
		// transform a 1D-array of real numbers
		// initializing
		int ndata = originalData.length;
		double[] alternatingOriginalData = new double[ndata*2];
		double[] alternatingTransformedData = new double[ndata*2];
		Complex[] transformedData = new Complex[ndata];

		// transforming
		for (int i = 0; i<ndata; i++) {
			alternatingOriginalData[2*i] = originalData[i];
			alternatingOriginalData[2*i+1] = 0.0D;
		}
		alternatingTransformedData = fft(alternatingOriginalData, 1);
		for (int i = 0; i<ndata; i++) {
			Complex ctemp = new Complex();
			ctemp.setReal(alternatingTransformedData[2*i]);
			ctemp.setImag(alternatingTransformedData[2*i+1]);
			transformedData[i] = ctemp;
		}

		return transformedData;
	}

	public static Complex[][] transform(double[][] originalData) {
		// transform a 2D-array of real numbers
		// initializing
		int m = originalData.length;
		int n = originalData[0].length;
		double[] alternatingOriginalData = new double[2*m*n];
		double[] alternatingTransformedData = new double[2*m*n];
		Complex[][] transformedData = new Complex[m][n];
		int[] arraySize = { m, n };

		// transforming
		int ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				alternatingOriginalData[ii] = originalData[j][i];
				alternatingOriginalData[ii+1] = 0.0D;
				ii += 2;
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 2, arraySize, 1);

		ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				Complex ctemp = new Complex();
				ctemp.setReal(alternatingTransformedData[ii]);
				ctemp.setImag(alternatingTransformedData[ii+1]);
				transformedData[j][i] = ctemp;
				ii += 2;
			}
		}

		return transformedData;
	}

	public static Complex[][][] transform(double[][][] originalData) {
		// transform a 3D-array of real numbers
		// initializing
		int l = originalData.length;
		int m = originalData[0].length;
		int n = originalData[0][0].length;
		double[] alternatingOriginalData = new double[2*m*n*l];
		double[] alternatingTransformedData = new double[2*m*n*l];
		Complex[][][] transformedData = new Complex[l][m][n];
		int[] arraySize = { l, m, n };

		// transforming
		int ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					alternatingOriginalData[ii] = originalData[k][j][i];
					alternatingOriginalData[ii+1] = 0.0D;
					ii += 2;
				}
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 3, arraySize, 1);

		ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					Complex ctemp = new Complex();
					ctemp.setReal(alternatingTransformedData[ii]);
					ctemp.setImag(alternatingTransformedData[ii+1]);
					transformedData[k][j][i] = ctemp;
					ii += 2;
				}
			}
		}

		return transformedData;
	}

	public static Complex[] transform(Complex[] originalData) {
		// transform a 1D-array of complex numbers
		// initializing
		int ndata = originalData.length;
		double[] alternatingOriginalData = new double[ndata*2];
		double[] alternatingTransformedData = new double[ndata*2];
		Complex[] transformedData = new Complex[ndata];

		// transforming
		for (int i = 0; i<ndata; i++) {
			alternatingOriginalData[2*i] = originalData[i].getReal();
			alternatingOriginalData[2*i+1] = originalData[i].getImag();
		}
		alternatingTransformedData = fft(alternatingOriginalData, 1);
		for (int i = 0; i<ndata; i++) {
			Complex ctemp = new Complex();
			ctemp.setReal(alternatingTransformedData[2*i]);
			ctemp.setImag(alternatingTransformedData[2*i+1]);
			transformedData[i] = ctemp;
		}

		return transformedData;
	}

	public static Complex[][] transform(Complex[][] originalData) {
		// transform a 2D-array of real numbers
		// initializing
		int m = originalData.length;
		int n = originalData[0].length;
		double[] alternatingOriginalData = new double[2*m*n];
		double[] alternatingTransformedData = new double[2*m*n];
		Complex[][] transformedData = new Complex[m][n];
		int[] arraySize = { m, n };

		// transforming
		int ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				alternatingOriginalData[ii] = originalData[j][i].getReal();
				alternatingOriginalData[ii+1] = originalData[j][i].getImag();
				ii += 2;
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 2, arraySize, 1);

		ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				Complex ctemp = new Complex();
				ctemp.setReal(alternatingTransformedData[ii]);
				ctemp.setImag(alternatingTransformedData[ii+1]);
				transformedData[j][i] = ctemp;
				ii += 2;
			}
		}

		return transformedData;
	}

	public static Complex[][][] transform(Complex[][][] originalData) {
		// transform a 3D-array of complex numbers
		// initializing
		int l = originalData.length;
		int m = originalData[0].length;
		int n = originalData[0][0].length;
		double[] alternatingOriginalData = new double[2*m*n*l];
		double[] alternatingTransformedData = new double[2*m*n*l];
		Complex[][][] transformedData = new Complex[l][m][n];
		int[] arraySize = { l, m, n };

		// transforming
		int ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					alternatingOriginalData[ii] = originalData[k][j][i].getReal();
					alternatingOriginalData[ii+1] = originalData[k][j][i].getImag();
					ii += 2;
				}
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 3, arraySize, 1);

		ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					Complex ctemp = new Complex();
					ctemp.setReal(alternatingTransformedData[ii]);
					ctemp.setImag(alternatingTransformedData[ii+1]);
					transformedData[k][j][i] = ctemp;
					ii += 2;
				}
			}
		}

		return transformedData;
	}

	public static Complex[] inverse(Complex[] originalData) {
		// inverse transform a 1D-array of Complex numbers
		// initializing
		int ndata = originalData.length;
		double[] alternatingOriginalData = new double[ndata*2];
		double[] alternatingTransformedData = new double[ndata*2];
		Complex[] transformedData = new Complex[ndata];

		// transforming
		for (int i = 0; i<ndata; i++) {
			alternatingOriginalData[2*i] = originalData[i].getReal();
			alternatingOriginalData[2*i+1] = originalData[i].getImag();
		}
		alternatingTransformedData = fft(alternatingOriginalData, -1);
		for (int i = 0; i<ndata; i++) {
			Complex ctemp = new Complex();
			ctemp.setReal(alternatingTransformedData[2*i]/ndata);
			ctemp.setImag(alternatingTransformedData[2*i+1]/ndata);
			transformedData[i] = ctemp;
		}

		return transformedData;
	}

	public static Complex[][] inverse(Complex[][] originalData) {
		// inverse transform a 2D-array of complex numbers
		// initializing
		int m = originalData.length;
		int n = originalData[0].length;
		double[] alternatingOriginalData = new double[2*m*n];
		double[] alternatingTransformedData = new double[2*m*n];
		Complex[][] transformedData = new Complex[m][n];
		int[] arraySize = { m, n };
		int ntot = m*n;

		// transforming
		int ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				alternatingOriginalData[ii] = originalData[j][i].getReal();
				alternatingOriginalData[ii+1] = originalData[j][i].getImag();
				ii += 2;
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 2, arraySize, -1);

		ii = 0;
		for (int j = 0; j<m; j++) {
			for (int i = 0; i<n; i++) {
				Complex ctemp = new Complex();
				ctemp.setReal(alternatingTransformedData[ii]/ntot);
				ctemp.setImag(alternatingTransformedData[ii+1]/ntot);
				transformedData[j][i] = ctemp;
				ii += 2;
			}
		}

		return transformedData;
	}

	public static Complex[][][] inverse(Complex[][][] originalData) {
		// inverse transform a 3D-array of complex numbers
		// initializing
		int l = originalData.length;
		int m = originalData[0].length;
		int n = originalData[0][0].length;
		double[] alternatingOriginalData = new double[2*m*n*l];
		double[] alternatingTransformedData = new double[2*m*n*l];
		Complex[][][] transformedData = new Complex[l][m][n];
		int[] arraySize = { l, m, n };
		int ntot = m*n*l;

		// transforming
		int ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					alternatingOriginalData[ii] = originalData[k][j][i].getReal();
					alternatingOriginalData[ii+1] = originalData[k][j][i].getImag();
					ii += 2;
				}
			}
		}

		alternatingTransformedData = fftn(alternatingOriginalData, 3, arraySize, -1);

		ii = 0;
		for (int k = 0; k<l; k++) {
			for (int j = 0; j<m; j++) {
				for (int i = 0; i<n; i++) {
					Complex ctemp = new Complex();
					ctemp.setReal(alternatingTransformedData[ii]/ntot);
					ctemp.setImag(alternatingTransformedData[ii+1]/ntot);
					transformedData[k][j][i] = ctemp;
					ii += 2;
				}
			}
		}

		return transformedData;
	}

	/**
     * Perform a Fast Fourier Transformation on 1D data array.
     * 
     * 
     * @param data data array to perform FFT on
     * @param isign 1 for regular FFT, -1 for inverse (result is NOT normalized
     * by N!)
     * 
     * @return string output
     */
	public static double[] fft(double[] data, int isign) {
		double dtemp = 0.0D, wtemp = 0.0D, tempr = 0.0D, tempi = 0.0D;
		double theta = 0.0D, wr = 0.0D, wpr = 0.0D, wpi = 0.0D, wi = 0.0D;
		int istep = 0, m = 0, mmax = 0;
		int nn = data.length/2;
		int n = nn<<1;
		int j = 1;

		for (int i = 1; i<n; i += 2) {
			if (j>i) {
				dtemp = data[j-1];
				data[j-1] = data[i-1];
				data[i-1] = dtemp;
				dtemp = data[j];
				data[j] = data[i];
				data[i] = dtemp;
			}
			m = n>>1;
			while (m>=2&&j>m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}
		mmax = 2;
		while (n>mmax) {
			istep = mmax<<1;
			theta = isign*(6.28318530717959D/mmax);
			wtemp = Math.sin(0.5D*theta);
			wpr = -2.0D*wtemp*wtemp;
			wpi = Math.sin(theta);
			wr = 1.0D;
			wi = 0.0D;
			for (m = 1; m<mmax; m += 2L) {
				for (int i = m; i<=n; i += istep) {
					j = i+mmax;
					tempr = wr*data[j-1]-wi*data[j];
					tempi = wr*data[j]+wi*data[j-1];
					data[j-1] = data[i-1]-tempr;
					data[j] = data[i]-tempi;
					data[i-1] += tempr;
					data[i] += tempi;
				}
				wr = (wtemp = wr)*wpr-wi*wpi+wr;
				wi = wi*wpr+wtemp*wpi+wi;
			}
			mmax = istep;
		}
		return data;
	}

	/**
     * Perform a Fast Fourier Transformation on n-D data array .
     * 
     * 
     * @param data data array to perform FFT on
     * @param dim number of dimensions
     * @param nn size of the array data in all directions
     * @param isign 1 for regular FFT, -1 for inverse (result is NOT normalized
     * by N!)
     * 
     * @return string output
     */
	public static double[] fftn(double[] data, int ndim, int[] nn, int isign) {
		int idim, i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
		int ibit, k1, k2, n, nprev, nrem, ntot;
		double tempi, tempr, dtemp, theta, wi, wpi, wpr, wr, wtemp;

		ntot = 1;
		for (idim = 0; idim<ndim; idim++) {
			ntot *= nn[idim];
		}
		nprev = 1;
		for (idim = ndim-1; idim>=0; idim--) {
			n = nn[idim];
			nrem = ntot/(n*nprev);
			ip1 = nprev<<1;
			ip2 = ip1*n;
			ip3 = ip2*nrem;
			i2rev = 0;
			for (i2 = 0; i2<ip2; i2 += ip1) {
				if (i2<=i2rev) {
					for (i1 = i2; i1<=i2+ip1-2; i1 += 2) {
						for (i3 = i1; i3<ip3; i3 += ip2) {
							i3rev = i2rev+i3-i2;
							dtemp = data[i3];
							data[i3] = data[i3rev];
							data[i3rev] = dtemp;
							dtemp = data[i3+1];
							data[i3+1] = data[i3rev+1];
							data[i3rev+1] = dtemp;
						}
					}
				}
				ibit = ip2>>1;
				while (ibit>=ip1&&i2rev>=ibit) {
					i2rev -= ibit;
					ibit >>= 1;
				}
				i2rev += ibit;
			}
			ifp1 = ip1;
			while (ifp1<ip2) {
				ifp2 = ifp1<<1;
				theta = isign*2*Math.PI/(ifp2/ip1);
				wtemp = Math.sin(0.5D*theta);
				wpr = -2.0D*wtemp*wtemp;
				wpi = Math.sin(theta);
				wr = 1.0D;
				wi = 0.0D;
				for (i3 = 0; i3<ifp1; i3 += ip1) {
					for (i1 = i3; i1<=i3+ip1-2; i1 += 2) {
						for (i2 = i1; i2<ip3; i2 += ifp2) {
							k1 = i2;
							k2 = k1+ifp1;
							tempr = wr*data[k2]-wi*data[k2+1];
							tempi = wr*data[k2+1]+wi*data[k2];
							data[k2] = data[k1]-tempr;
							data[k2+1] = data[k1+1]-tempi;
							data[k1] += tempr;
							data[k1+1] += tempi;
						}
					}
					wr = (wtemp = wr)*wpr-wi*wpi+wr;
					wi = wi*wpr+wtemp*wpi+wi;
				}
				ifp1 = ifp2;
			}
			nprev *= n;
		}
		return data;
	}

	/**
     * checks whether n is an integer power of 2
     * @param n
     * @return
     */
	public static boolean checkPowerOfTwo(int n) {
		int temp = n;
		while (temp>1) {
			if ((temp%2)!=0) {
				return false;
			} else {
				temp /= 2;
			}
		}
		return true;
	}

}