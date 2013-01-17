/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package simulator.reaction.kinetic;

import org.jdom.Element;

import utils.XMLParser;

public class LinearKinetic extends IsKineticFactor {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private double _K;
	
	public LinearKinetic() {
	}

	public void init(Element defMarkUp) {
		_K = (new XMLParser(defMarkUp)).getParamDbl("K");
		nParam = 1;
	}

	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex) {
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("K");
	}

	public double kineticValue(double solute, double[] paramTable, int index) {
		return paramTable[index]*solute;
	}

	public double kineticValue(double solute) {
		return _K*solute;
	}

	public double kineticDiff(double solute, double[] paramTable, int index) {
		return paramTable[index];
	}

	public double kineticDiff(double solute) {
		return _K;
	}

	public double kineticMax() {
		return 1;
	}

	public double kineticMax(double[] paramTable, int index) {
		return 1;
	}
}
