/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *______________________________________________________
 * Implements a reference to a grid element
 * 
 */

/**
 * @since January 2004
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 */

package simulator.geometry;

public class GridElementFloat {

	/**
     * Sorts the grid elements according to the T value
     */
	public static class TValueComparator implements java.util.Comparator<Object> {

		public int compare(Object b1, Object b2) {
			double z1 = ((GridElementFloat) b1).getValue();
			double z2 = ((GridElementFloat) b2).getValue();
			return (z1>z2 ? 1 : -1);
		}
	}

	public int _index;
	private double[] _values;

	/**
     * Initialize a new grid element
     * 
     * @param i
     * @param j
     * @param k
     * @param v
     */
	public GridElementFloat(double[] values, int index) {
		_index = index;
		_values = values;
	}

	/**
     * @return the levelset value for this grid element
     */
	public double getValue() {
		return _values[_index];
	}
}