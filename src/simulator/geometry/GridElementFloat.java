/**
 * \package simulator.geometry
 * \brief Package of boundary utilities that aid the creation of the environment being simulated
 * 
 * Package of boundary utilities that aid the creation of the environment being simulated. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry;

/**
 * \brief Implements a reference to a grid element
 * 
 * Implements a reference to a grid element
 * 
 * @author Joao Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public class GridElementFloat 
{
	/**
	 * Index to a location on a grid which relates to this element
	 */
	public int _index;
	
	/**
	 * Array of double values that comprise this grid
	 */
	private double[] _values;

	/**
     * Initialize a new grid element
     * 
     * @param values	Array of doubles to set the grid to
     * @param index	Integer index to set the index to
     */
	public GridElementFloat(double[] values, int index) {
		_index = index;
		_values = values;
	}

	/**
     * \brief Return the double value in the array of doubles for the specified index
     * 
     * Return the double value in the array of doubles for the specified index
     * 
     * @return Double value stored in the array
     */
	public double getValue() {
		return _values[_index];
	}
	
	/**
     * \brief Sorts the grid elements according to the T value
     * 
     * Sorts the grid elements according to the T value
     */
	public static class TValueComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			double z1 = ((GridElementFloat) b1).getValue();
			double z2 = ((GridElementFloat) b2).getValue();
			return (z1>z2 ? 1 : -1);
		}
	}
}