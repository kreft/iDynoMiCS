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
 * \brief Class for storing an 1- to 3-dimensional array of complex numbers as an 1D array of alternating real and imaginary parts in double precision. 
 * 
 * Class for storing an 1- to 3-dimensional array of complex numbers as an 1D array of alternating real and imaginary parts in double precision. 
 * This data structure optimizes efficiency of fast fourier transformation (fft) by the "Numerical Recipes" algorithm. Includes basic 
 * math operations (addition, subtraction, multiplication, division)
 * 
 * @author Andreas D�tsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 *
 */
public class ComplexArray {

	/**
	 * The main array
	 */
	protected double[]	altArray;
	
	/**
	 * Dimension of the source array
	 */
	protected int		_dim;
	
	/**
	 * Number of voxels in the l direction
	 */
	protected int		_l;
	
	/**
	 * Number of voxels in the m direction
	 */
	protected int 		_m;
	
	/**
	 * Number of voxels in the n direction
	 */
	protected int		_n;
	
	/**
	 * Integer array storing the original size of the array (nxmxl)
	 */
	protected int[]		_originalSize	= {_n, _m, _l};
	
	/**
	 *  Length of the main array
	 */
	protected int		_ntot;

	/**
	 * \brief Constructor to generate a 1D array of complex numbers
	 * 
	 * Constructor to generate a 1D array of complex numbers
	 * 
	 * @param n	The number of entries that will be included in this array
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
	 * \brief Constructor to generate a 2D array of complex numbers
	 * 
	 * Constructor to generate a 2D array of complex numbers
	 * 
	 * @param n	The number of entries in the n direction that will be included in this array
	 * @param m	The number of entries in the m direction that will be included in this array
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
	 * \brief Constructor to generate a 3D array of complex numbers
	 * 
	 * Constructor to generate a 3D array of complex numbers
	 * 
	 * @param n	The number of entries in the n direction that will be included in this array
	 * @param m	The number of entries in the m direction that will be included in this array
	 * @param L	The number of entries in the L direction that will be included in this array
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

	/**
	 * \brief Return the dimension of this complex array
	 * 
	 * Return the dimension of this complex array
	 * 
	 * @return	Dimension of this array
	 */
	public int getDimension() {
		return this._dim;
	}

	/**
	 * \brief Return the total length of this complex array
	 * 
	 * Return the total length of this complex array
	 * 
	 * @return	Total length of this array
	 */
	public int getTotalLength() {
		return this._ntot;
	}

	/**
	 * \brief Return the array that shows the size of the original array before storing using this method
	 * 
	 * Return the array that shows the size of the original array before storing using this method
	 * 
	 * @return	Array showing the original nxmxl size of the array
	 */
	public int[] getOriginalSize() {
		return this._originalSize;
	}

	/**
	 * \brief Return the complex array
	 * 
	 * Return the complex array
	 * 
	 * @return Complex array (of doubles)
	 */
	public double[] getWholeArray() {
		return this.altArray;
	}

	/**
	 * \brief Set the array to values of a array of doubles
	 * 
	 * Set the array to values of a array of doubles
	 * 
	 * @param inputArray	Array of values to set this complex array to
	 */
	public void setWholeArray(double[] inputArray) {
		for (int i = 0; i < inputArray.length; i++) {
			this.altArray[i] = inputArray[i];
		}
	}

	/**
	 * \brief Set the array to values to those of another AltComplex
	 * 
	 * Set the array to values to those of another AltComplex
	 * 
	 * @param input	Complex array of values whose values should be cloned into this array
	 */
	public void setWholeArray(ComplexArray input) {
		double[] temp = input.getWholeArray();
		for (int i = 0; i < temp.length; i++) {
			this.altArray[i] = temp[i];
		}
	}

	/**
	 * \brief Set value at given coordinates
	 * 
	 * Set value at given coordinates
	 * 
	 * @param realPart	The value to set the array coordinate to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void setRealAt(double realPart, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] = realPart;
	}

	/**
	 * \brief Set value at given coordinates
	 * 
	 * Set value at given coordinates
	 * 
	 * @param realPart	The value to set the array coordinate to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void setRealAt(double realPart, int i, int j) {
		this.setRealAt(realPart, i, j, 0);
	}
	
	/**
	 * \brief Set value at given coordinates
	 * 
	 * Set value at given coordinates
	 * 
	 * @param realPart	The value to set the array coordinate to
	 * @param i	Array coordinate i
	 */ 
	public void setRealAt(double realPart, int i) {
		this.setRealAt(realPart, i, 0, 0);
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a value
	 * 
	 * Set the 'imaginary' part of the array to a value
	 * 
	 * @param imagPart	Value to set this part of the array to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void setImagAt(double imagPart, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n) + 1;
		this.altArray[q] = imagPart;
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a value
	 * 
	 * Set the 'imaginary' part of the array to a value
	 * 
	 * @param imagPart	Value to set this part of the array to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void setImagAt(double imagPart, int i, int j) {
		this.setImagAt(imagPart, i, j, 0);
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a value
	 * 
	 * Set the 'imaginary' part of the array to a value
	 * 
	 * @param imagPart	Value to set this part of the array to
	 * @param i	Array coordinate i
	 */
	public void setImagAt(double imagPart, int i) {
		this.setImagAt(imagPart, i, 0, 0);
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a Complex number
	 * 
	 * Set the 'imaginary' part of the array to a Complex number
	 * 
	 * @param c	Complex number to set this part of the array to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void setComplexAt(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] = c.getReal();
		this.altArray[q + 1] = c.getImag();
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a Complex number
	 * 
	 * Set the 'imaginary' part of the array to a Complex number
	 * 
	 * @param c	Complex number to set this part of the array to
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void setComplexAt(Complex c, int i, int j) {
		this.setComplexAt(c, i, j, 0);
	}

	/**
	 * \brief Set the 'imaginary' part of the array to a Complex number
	 * 
	 * Set the 'imaginary' part of the array to a Complex number
	 * 
	 * @param c	Complex number to set this part of the array to
	 * @param i	Array coordinate i
	 */
	public void setComplexAt(Complex c, int i) {
		this.setComplexAt(c, i, 0, 0);
	}

	/**
	 * \brief Get value at specified coordinates
	 * 
	 * Get value at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 * @return	Value at that coordinate
	 */
	public double getRealAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		return this.altArray[q];
	}

	/**
	 * \brief Get value at specified coordinates
	 * 
	 * Get value at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @return	Value at that coordinate
	 */
	public double getRealAt(int i, int j) {
		int q = 2 * (i + j * _n);
		return this.altArray[q];
	}

	/**
	 * \brief Get value at specified coordinates
	 * 
	 * Get value at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @return	Value at that coordinate
	 */
	public double getRealAt(int i){
		return this.altArray[2 * i];
	}

	/**
	 * \brief Get imaginary value part of this object at specified coordinates
	 * 
	 * Get imaginary value part of this object at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 * @return	Value at that coordinate
	 */
	public double getImagAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n) + 1;
		return this.altArray[q];
	}

	/**
	 * \brief Get imaginary value part of this object at specified coordinates
	 * 
	 * Get imaginary value part of this object at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @return	Value at that coordinate
	 */
	public double getImagAt(int i, int j) {
		int q = 2 * (i + j * _n) + 1;
		return this.altArray[q];
	}

	/**
	 * \brief Get imaginary value part of this object at specified coordinates
	 * 
	 * Get imaginary value part of this object at specified coordinates
	 * 
	 * @param i	Array coordinate i
	 * @return	Value at that coordinate
	 */
	public double getImagAt(int i) {
		return this.altArray[2 * i + 1];
	}

	/**
	 * \brief Get the complex number stored at a given position in the array
	 * 
	 * Get the complex number stored at a given position in the array
	 * 
	 * @param i	Array coordinate i
	 * @param j Array coordinate j
	 * @param k	Array coordinate k
	 * @return	Complex number at that position
	 */
	public Complex getComplexAt(int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		return new Complex(this.altArray[q], this.altArray[q + 1]);
	}

	/**
	 * \brief Get the complex number stored at a given position in the array
	 * 
	 * Get the complex number stored at a given position in the array
	 * 
	 * @param i	Array coordinate i
	 * @param j Array coordinate j
	 * @return	Complex number at that position
	 */
	public Complex getComplexAt(int i, int j) {
		int q = 2 * (i + j * _n);
		return new Complex(this.altArray[q], this.altArray[q + 1]);
	}

	/**
	 * \brief Get the complex number stored at a given position in the array
	 * 
	 * Get the complex number stored at a given position in the array
	 * 
	 * @param i	Array coordinate i
	 * @return	Complex number at that position
	 */
	public Complex getComplexAt(int i) {
		return new Complex(this.altArray[2 * i], this.altArray[2 * i + 1]);
	}

	/**
	 * \brief Add a double to the data in a given position in the array
	 * 
	 * Add a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void add(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] += realNumber;
	}

	/**
	 * \brief Subtract a double to the data in a given position in the array
	 * 
	 * Subtract a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to subtract
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void sub(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] -= realNumber;
	}

	/**
	 * \brief Multiply a double to the data in a given position in the array
	 * 
	 * Multiply a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to multiply
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void mul(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] *= realNumber;
		this.altArray[q + 1] *= realNumber;
	}

	/**
	 * \brief Divide a double to the data in a given position in the array
	 * 
	 * Divide a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to divide by
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void div(double realNumber, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] /= realNumber;
		this.altArray[q + 1] /= realNumber;
	}

	/**
	 * \brief Add a double to the data in a given position in the array
	 * 
	 * Add a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void add(double realNumber, int i, int j) {
		this.add(realNumber, i, j, 0);
	}

	/**
	 * \brief Subtract a double to the data in a given position in the array
	 * 
	 * Subtract a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to subtract
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void sub(double realNumber, int i, int j) {
		this.sub(realNumber, i, j, 0);
	}

	/**
	 * \brief Multiply a double to the data in a given position in the array
	 * 
	 * Multiply a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to multiply
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void mul(double realNumber, int i, int j) {
		this.mul(realNumber, i, j, 0);
	}

	/**
	 * \brief Divide a double to the data in a given position in the array
	 * 
	 * Divide a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to divide by
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void div(double realNumber, int i, int j) {
		this.div(realNumber, i, j, 0);
	}

	/**
	 * \brief Add a double to the data in a given position in the array
	 * 
	 * Add a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to add
	 * @param i	Array coordinate i
	 */
	public void add(double realNumber, int i) {
		this.add(realNumber, i, 0, 0);
	}

	/**
	 * \brief Subtract a double to the data in a given position in the array
	 * 
	 * Subtract a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to subtract
	 * @param i	Array coordinate i
	 */
	public void sub(double realNumber, int i) {
		this.sub(realNumber, i, 0, 0);
	}

	/**
	 * \brief Multiply a double to the data in a given position in the array
	 * 
	 * Multiply a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to multiply
	 * @param i	Array coordinate i
	 */
	public void mul(double realNumber, int i) {
		this.mul(realNumber, i, 0, 0);
	}

	/**
	 * \brief Divide a double to the data in a given position in the array
	 * 
	 * Divide a double to the data in a given position in the array
	 * 
	 * @param realNumber Value to divide by
	 * @param i	Array coordinate i
	 */
	public void div(double realNumber, int i) {
		this.div(realNumber, i, 0, 0);
	}

	/**
	 * \brief Add a complex number to the data in a given position in the array
	 * 
	 * Add a complex number to the data in a given position in the array
	 * 
	 * @param c Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void add(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] += c.getReal();
		this.altArray[q + 1] += c.getImag();
	}

	/**
	 * \brief Subtract a complex number from the data in a given position in the array
	 * 
	 * Subtract a complex number from the data in a given position in the array
	 * 
	 * @param c Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void sub(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		this.altArray[q] -= c.getReal();
		this.altArray[q + 1] -= c.getImag();
	}

	/**
	 * \brief Multiply a complex number to the data in a given position in the array
	 * 
	 * Multiply a complex number to the data in a given position in the array
	 * 
	 * @param c Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
	public void mul(Complex c, int i, int j, int k) {
		int q = 2 * (i + j * _n + k * _m * _n);
		double cRe = c.getReal();
		double cIm = c.getImag();
		double temp = this.altArray[q] * cRe - this.altArray[q + 1] * cIm;
		this.altArray[q + 1] = this.altArray[q + 1] * cRe + this.altArray[q]
				* cIm;
		this.altArray[q] = temp;
	}

	/**
	 * \brief Divide a complex number to the data in a given position in the array
	 * 
	 * Divide a complex number to the data in a given position in the array
	 * 
	 * @param c Value to divide by
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 * @param k	Array coordinate k
	 */
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

	/**
	 * \brief Add a complex number to the data in a given position in the array
	 * 
	 * Add a complex number to the data in a given position in the array
	 * 
	 * @param c Value to add
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void add(Complex c, int i, int j) {
		this.add(c, i, j, 0);
	}

	/**
	 * \brief Add a complex number to the data in a given position in the array
	 * 
	 * Add a complex number to the data in a given position in the array
	 * 
	 * @param c Value to add
	 * @param i	Array coordinate i
	 */
	public void add(Complex c, int i) {
		this.add(c, i, 0, 0);
	}


	/**
	 * \brief Subtract a complex number from the data in a given position in the array
	 * 
	 * Subtract a complex number from the data in a given position in the array
	 * 
	 * @param c Value to subtract
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void sub(Complex c, int i, int j) {
		this.sub(c, i, j, 0);
	}

	/**
	 * \brief Subtract a complex number from the data in a given position in the array
	 * 
	 * Subtract a complex number from the data in a given position in the array
	 * 
	 * @param c Value to subtract
	 * @param i	Array coordinate i
	 */
	public void sub(Complex c, int i) {
		this.sub(c, i, 0, 0);
	}

	/**
	 * \brief Multiply a complex number to the data in a given position in the array
	 * 
	 * Multiply a complex number to the data in a given position in the array
	 * 
	 * @param c Value to multiply
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void mul(Complex c, int i, int j) {
		this.mul(c, i, j, 0);
	}

	/**
	 * \brief Multiply a complex number to the data in a given position in the array
	 * 
	 * Multiply a complex number to the data in a given position in the array
	 * 
	 * @param c Value to multiply
	 * @param i	Array coordinate i
	 */
	public void mul(Complex c, int i) {
		this.mul(c, i, 0, 0);
	}

	/**
	 * \brief Divide a complex number to the data in a given position in the array
	 * 
	 * Divide a complex number to the data in a given position in the array
	 * 
	 * @param c Value to divide
	 * @param i	Array coordinate i
	 * @param j	Array coordinate j
	 */
	public void div(Complex c, int i, int j) {
		this.div(c, i, j, 0);
	}

	/**
	 * \brief Divide a complex number to the data in a given position in the array
	 * 
	 * Divide a complex number to the data in a given position in the array
	 * 
	 * @param c Value to divide by
	 * @param i	Array coordinate i
	 */
	public void div(Complex c, int i) {
		this.div(c, i, 0, 0);
	}

}
