/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

package simulator.reaction.kinetic;

import org.jdom.Element;

import utils.ExtraMath;
import utils.XMLParser;

public class MonodKinetic extends IsKineticFactor {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private double _Ks;

	// Not called
	public MonodKinetic() {
	}

	// Not called
	public MonodKinetic(double Ks) {
		_Ks = Ks;
		nParam = 1;
	}

	public void init(Element defMarkUp) {
		_Ks = (new XMLParser(defMarkUp)).getParamDbl("Ks");
		nParam = 1;
	}

	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex) {
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("Ks");
	}

	public double kineticValue(double solute, double[] paramTable, int index) {
		return solute/(paramTable[index]+solute);
	}

	public double kineticValue(double solute) {
		return solute/(_Ks+solute);
	}

	public double kineticDiff(double solute, double[] paramTable, int index) {
		return paramTable[index]/ExtraMath.sq(paramTable[index]+solute);
	}

	public double kineticDiff(double solute) {
		return _Ks/ExtraMath.sq(_Ks+solute);
	}

	// Not called
	public double kineticMax() {
		return 1;
	}

	// Not called
	public double kineticMax(double[] paramTable, int index) {
		return 1;
	}
}
