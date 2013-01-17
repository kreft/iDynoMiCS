/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package simulator.reaction.kinetic;

import org.jdom.Element;

/**
 * Rob 7/12/2011
 * Use of this kinetic is discouraged: as the kinetic value is independent of 
 * solute concentration this can lead to negative concentration values. Note 
 * that LinearKinetic is a more suitable choice in many situations.
 */
public class FirstOrderKinetic extends IsKineticFactor {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	public FirstOrderKinetic() {
	}

	public void init(Element defMarkUp) {
		nParam = 0;
	}

	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex) {
	}

	public double kineticValue(double solute, double[] paramTable, int index) {
		return 1;
	}

	public double kineticValue(double solute) {
		return 1;
	}

	public double kineticDiff(double solute, double[] paramTable, int index) {
		return 0;
	}

	public double kineticDiff(double solute) {
		return 0;
	}

	public double kineticMax() {
		return 1;
	}

	public double kineticMax(double[] paramTable, int index) {
		return 1;
	}
}
