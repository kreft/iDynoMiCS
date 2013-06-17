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
 * \brief Creates a 2D or 3D Vectorized array of a set size, and provides utilities to set and get values from this array
 * 
 * Creates a 2D or 3D Vectorized array of a set size, and provides utilities to set and get values from this array
 */
public class VectorizedArray 
{
	/**
	 * Vectorized array to create  
	 */
	private double[] _array;

	/**
	 * I (X) size of original array
	 */
	private int _nI;
	
	/**
	 * J (J) size of original array
	 */
	private int _nJ;
	
	/**
	 * K (Z) size of original array
	 */
	private int _nK;

	// 
	@SuppressWarnings("unused")
	/**
	 * Dimensionality of original array
	 */
	private int _dim;

	/**
	 * \brief Creates an empty 3D array of three provided dimensions
	 * 
	 * Creates an empty 3D array of three provided dimensions
	 * 
	 * @param nI	Array size in I (X) direction
	 * @param nJ	Array size in J (Y) direction
	 * @param nK	Array size in K (Z) direction
	 */
	public VectorizedArray(int nI, int nJ, int nK) {
		_nI = nI;
		_nJ = nJ;
		_nK = nK;
		_dim = 3;
		_array = new double[nI * nJ * nK];
	}

	/**
	 * \brief Creates an empty 2D array in two provided dimensions
	 * 
	 * Creates an empty 2D array in two provided dimensions
	 * 
	 * @param nI	Array size in I (X) direction
	 * @param nJ	Array size in J (Y) direction
	 */
	public VectorizedArray(int nI, int nJ) {
		_nI = nI;
		_nJ = nJ;
		_nK = 1;
		_dim = 2;
		_array = new double[nI * nJ];
	}

	/**
	 * \brief Return the value at a specified index in the array
	 * 
	 * Return the value at a specified index in the array
	 * 
	 * @param i	Integer index of the i position for which the value is being sought
	 * @param j Integer index of the j position for which the value is being sought
	 * @param k Integer index of the k position for which the value is being sought
	 * @return 	Double value at this position
	 */
	public double getValueAt(int i, int j, int k) {
		return _array[i + _nI * j + _nI * _nJ * k];
	}

	/**
	 * \brief Set a value at a given array position
	 * 
	 * Set a value at a given array position
	 * 
	 * @param value	Value to set position to
	 * @param i	Integer index of the i position for which the value is to be set
	 * @param j Integer index of the j position for which the value is to be set
	 * @param k Integer index of the k position for which the value is to be set     
	 */
	public void setValueAt(double value, int i, int j, int k) {
		_array[i + _nI * j + _nI * _nJ * k] = value;
	}

	/**
	 * \brief Return the maximum value of the array
	 * 
	 * Return the maximum value of the array
	 * 
	 * @return Double value which is the maximum value in the array
	 */
	public double getMax() {
		double maxValue = _array[0];
		for (int i = 1; i < _array.length; i++) {
			if (_array[i] > maxValue)
				maxValue = _array[i];
		}
		return maxValue;
	}
	
	/**
	 * \brief Return the minimum value of the array
	 * 
	 * Return the minimum value of the array
	 * 
	 * @return Double value which is the minimum value in the array
	 */
	public double getMin() {
		double minValue = _array[0];
		for (int i = 1; i < _array.length; i++) {
			if (_array[i] < minValue)
				minValue = _array[i];
		}
		return minValue;
	}

	/**
	 * \brief Return the array as a 3D double array
	 * 
	 * Return the array as a 3D double array
	 * 
	 * @return Vectorized array as a 3D double array
	 */
	public double[][][] getValues() {
		double[][][] array = new double[_nI][_nJ][_nK];
		for (int i = 0; i < _nI; i++) {
			for (int j = 0; j < _nJ; j++) {
				for (int k = 0; k < _nK; k++) {
					array[i][j][k] = _array[i + _nI * j + _nI * _nJ * k];
				}
			}
		}
		return array;
	}
}
