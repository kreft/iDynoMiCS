/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 * ______________________________________________________
 * Class for storing an 1- to 3-dimensional array of complex numbers as an 1D
 * array of alternating real and imaginary parts in double precision. This data
 * structure optimizes efficiency of fast fourier transformation (fft) by the
 * "Numerical Recipes" algorithm. Includes basic math operations (addition,
 * subtraction, multiplication, division)
 * 
 */

/**
 * 
 * @since August 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 */
 
package utils;
 
public class ComplexArray {

	// the main array
	protected double[]	altArray;
	// dimension of the source array
	protected int		_dim;
	// number of voxels in each direction
	protected int		_l, _m, _n;
	protected int[]		_originalSize	= {_n, _m, _l};
	// length of the altArray
	protected int		_ntot;

	/**
	 * constructor for an 1D array of Complex numbers
	 * 
	 * @param n
	 */
	public ComplexArray(int n) {
		_n = n;
		_m = 1;
		_l = 1;
		_dim = 1;
		_originalSize[0] = _n;
		_originalSize[1] = _m;
		_originalSize[2] = _l;
		altArray = new double[2 * _n];
	}

	/**
	 * constructor for an 2D array of Complex numbers
	 * 
	 * @param n
	 */
	public ComplexArray(int n, int m) {
		_n = n;
		_m = m;
		_l = 1;
		_dim = 2;
		_originalSize[0] = _n;
		_originalSize[1] = _m;
		_originalSize[2] = _l;
		altArray = new double[2 * _n * _m];
	}

	/**
	 * constructor for an 3D array of Complex numbers
	 * 
	 * @param n
	 */
	public ComplexArray(int n, int m, int L) {
		_n = n;
		_m = m;
		_l = L;
		_dim = 3;
		_originalSize[0] = _n;
		_originalSize[1] = _m;
		_originalSize[2] = _l;
		altArray = new double[2 * _n * _m * _l];
	}

	public int getDimension() {
		return this._dim;
	}

	public int getTotalLength() {
		return this._ntot;
	}

	public int[] getOriginalSize() {
		return this._originalSize;
	}

	public double[] getWholeArray() {
		return this.altArray;
	}

	/**
	 * set the array to values of a array of doubles
	 * 
	 * @param inputArray
	 */
	public void setWholeArray(double[] inputArray) {
		for (int i = 0; i < inputArray.length; i++) {
			this.altArray[i] = inputArray[i];
		}
	}

	/**
	 * the array to values of another AltComplex
	 * 
	 * @param input
	 */
	public void setWholeArray(ComplexArray input) {
		double[] temp = input.getWholeArray();
		for (int i = 0; i < temp.length; i++) {
			this.altArray[i] = temp[i];
		}
	}

	// setting values at coordinates
	public void setRealAt(double realPart, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] = realPart;
	}

	public void setRealAt(double realPart, int i, int j) {
		this.setRealAt(realPart, i, j, 0);
	}

	public void setRealAt(double realPart, int i) {
		this.setRealAt(realPart, i, 0, 0);
	}

	public void setImagAt(double imagPart, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n) + 1;
		this.altArray[q] = imagPart;
	}

	public void setImagAt(double imagPart, int i, int j) {
		this.setImagAt(imagPart, i, j, 0);
	}

	public void setImagAt(double imagPart, int i) {
		this.setImagAt(imagPart, i, 0, 0);
	}

	public void setComplexAt(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] = c.getReal();
		this.altArray[q + 1] = c.getImag();
	}

	public void setComplexAt(Complex c, int i, int j) {
		this.setComplexAt(c, i, j, 0);
	}

	public void setComplexAt(Complex c, int i) {
		this.setComplexAt(c, i, 0, 0);
	}

	// getting values at coordinates
	public double getRealAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		return this.altArray[q];
	}

	public double getRealAt(int i, int j) {
		int q = 2 * (i + j * _n);
		return this.altArray[q];
	}

	public double getRealAt(int i) {
		return this.altArray[2 * i];
	}

	public double getImagAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n) + 1;
		return this.altArray[q];
	}

	public double getImagAt(int i, int j) {
		int q = 2 * (i + j * _n) + 1;
		return this.altArray[q];
	}

	public double getImagAt(int i) {
		return this.altArray[2 * i + 1];
	}

	public Complex getComplexAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		return new Complex(this.altArray[q], this.altArray[q + 1]);
	}

	public Complex getComplexAt(int i, int j) {
		int q = 2 * (i + j * _n);
		return new Complex(this.altArray[q], this.altArray[q + 1]);
	}

	public Complex getComplexAt(int i) {
		return new Complex(this.altArray[2 * i], this.altArray[2 * i + 1]);
	}

	// ### operators for doubles
	public void add(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] += realNumber;
	}

	public void sub(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] -= realNumber;
	}

	public void mul(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] *= realNumber;
		this.altArray[q + 1] *= realNumber;
	}

	public void div(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] /= realNumber;
		this.altArray[q + 1] /= realNumber;
	}

	public void add(double realNumber, int i, int j) {
		this.add(realNumber, i, j, 0);
	}

	public void sub(double realNumber, int i, int j) {
		this.sub(realNumber, i, j, 0);
	}

	public void mul(double realNumber, int i, int j) {
		this.mul(realNumber, i, j, 0);
	}

	public void div(double realNumber, int i, int j) {
		this.div(realNumber, i, j, 0);
	}

	public void add(double realNumber, int i) {
		this.add(realNumber, i, 0, 0);
	}

	public void sub(double realNumber, int i) {
		this.sub(realNumber, i, 0, 0);
	}

	public void mul(double realNumber, int i) {
		this.mul(realNumber, i, 0, 0);
	}

	public void div(double realNumber, int i) {
		this.div(realNumber, i, 0, 0);
	}

	// ### operators for complex numbers
	public void add(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] += c.getReal();
		this.altArray[q + 1] += c.getImag();
	}

	public void sub(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] -= c.getReal();
		this.altArray[q + 1] -= c.getImag();
	}

	public void mul(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		double cRe = c.getReal();
		double cIm = c.getImag();
		double temp = this.altArray[q] * cRe - this.altArray[q + 1] * cIm;
		this.altArray[q + 1] = this.altArray[q + 1] * cRe + this.altArray[q]
				* cIm;
		this.altArray[q] = temp;
	}

	public void div(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		double cRe = c.getReal();
		double cIm = c.getImag();
		double den = cRe * cRe + cIm * cIm;
		double temp = (this.altArray[q] * cRe + this.altArray[q + 1] * cIm)
				/ den;
		this.altArray[q + 1] = (this.altArray[q + 1] * cRe - this.altArray[q]
				* cIm)
				/ den;
		this.altArray[q] = temp;
	}

	public void add(Complex c, int i, int j) {
		this.add(c, i, j, 0);
	}

	public void add(Complex c, int i) {
		this.add(c, i, 0, 0);
	}

	public void sub(Complex c, int i, int j) {
		this.sub(c, i, j, 0);
	}

	public void sub(Complex c, int i) {
		this.sub(c, i, 0, 0);
	}

	public void mul(Complex c, int i, int j) {
		this.mul(c, i, j, 0);
	}

	public void mul(Complex c, int i) {
		this.mul(c, i, 0, 0);
	}

	public void div(Complex c, int i, int j) {
		this.div(c, i, j, 0);
	}

	public void div(Complex c, int i) {
		this.div(c, i, 0, 0);
	}

}
