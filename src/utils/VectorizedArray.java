
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *
 */

/**
 * @since June 2008
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 */

package utils;

public class VectorizedArray {

	// values
	private double[] _array;

	// size of original array
	private int _nI, _nJ, _nK;

	// dimensionality of orignal array
	@SuppressWarnings("unused")
	private int _dim;

	/**
	 * creates an empty 3D array
	 * 
	 * @param nI,nJ,nK
	 *            the size of the original array
	 */
	public VectorizedArray(int nI, int nJ, int nK) {
		_nI = nI;
		_nJ = nJ;
		_nK = nK;
		_dim = 3;
		_array = new double[nI * nJ * nK];
	}

	/**
	 * creates an empty 2D array
	 * 
	 * @param nI,nJ
	 *            the size of the original array
	 */
	public VectorizedArray(int nI, int nJ) {
		_nI = nI;
		_nJ = nJ;
		_nK = 1;
		_dim = 2;
		_array = new double[nI * nJ];
	}

	/**
	 * @param i,j,k
	 *            the original index
	 * @return the value at this position
	 */
	public double getValueAt(int i, int j, int k) {
		return _array[i + _nI * j + _nI * _nJ * k];
	}

	/**
	 * set a value at a certain position
	 * 
	 * @param i,j,k
	 *            the original index
	 */
	public void setValueAt(double value, int i, int j, int k) {
		_array[i + _nI * j + _nI * _nJ * k] = value;
	}

	/**
	 * @return the maximum value of the array
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
	 * @return the minimum value of the array
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
	 * @return the array as 3D double array
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
