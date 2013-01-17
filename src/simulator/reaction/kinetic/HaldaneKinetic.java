
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

package simulator.reaction.kinetic;

import org.jdom.Element;

import utils.ExtraMath;
import utils.XMLParser;
public class HaldaneKinetic extends IsKineticFactor {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private double _Ks;
	private double _Ki;

	public HaldaneKinetic(double Ks, double Ki) {
		_Ks = Ks;
		_Ki = Ki;
		nParam = 2;
	}

	public void init(Element defMarkUp) {
		_Ks = (new XMLParser(defMarkUp)).getParamDbl( "Ks");
		_Ki = (new XMLParser(defMarkUp)).getParamDbl( "Ki");
		nParam = 2;
	}

	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex) {
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("Ks");
		kineticParam[paramIndex+1] = (new XMLParser(defMarkUp)).getParamDbl("Ki");
	}

	public double kineticValue(double solute, double[] paramTable, int index) {
		return solute/(paramTable[index]+solute+solute*solute/paramTable[index+1]);
	}

	public double kineticValue(double solute) {
		return solute/(_Ks+solute+solute*solute/_Ki);
	}

	public double kineticDiff(double solute) {
		return (_Ks-ExtraMath.sq(solute)/_Ki)/ExtraMath.sq(_Ks+solute+solute*solute/_Ki);
	}

	public double kineticDiff(double solute, double[] paramTable, int index) {
		return (paramTable[index]-ExtraMath.sq(solute)/paramTable[index+1])
		        /ExtraMath.sq(paramTable[index]+solute+solute*solute/paramTable[index+1]);
	}

	public double kineticMax() {
		double S_null = Math.sqrt(_Ks*_Ki);
		return kineticValue(S_null);
	}

	public double kineticMax(double[] paramTable, int index) {
		double S_null = Math.sqrt(_Ks*_Ki);
		return kineticValue(S_null, paramTable, index);
	}
}
