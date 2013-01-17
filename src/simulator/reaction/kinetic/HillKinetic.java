/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

package simulator.reaction.kinetic;

import org.jdom.Element;
import utils.ExtraMath;
import utils.XMLParser;


public class HillKinetic extends IsKineticFactor {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	// parameters
	private double            _Ks;
	private double            _h;
	// auxiliaries
	private double            _KsH, _KsPowH;

	public HillKinetic() {
	}
	
	public HillKinetic(double Ks, double h) {
		_Ks = Ks;
		_h = h;
		_KsH = _Ks*_h;
		_KsPowH = Math.pow(_Ks, _h);
		nParam = 2;
	}

	public void init(Element defMarkUp) {
		_Ks = (new XMLParser(defMarkUp)).getParamDbl("Ks");
		_h = (new XMLParser(defMarkUp)).getParamDbl("h");
		_KsH = Math.pow(_Ks, _h)*_h;
		_KsPowH = Math.pow(_Ks, _h);
		nParam = 2;
	}

	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex) {
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl( "Ks");
		kineticParam[paramIndex+1] = (new XMLParser(defMarkUp)).getParamDbl("h");
	}

	public double kineticValue(double solute, double[] paramTable, int index) {
		return Math.pow(solute, paramTable[index+1])
		        /(Math.pow(paramTable[index], paramTable[index+1])+Math.pow(solute,
		                paramTable[index+1]));
	}

	public double kineticValue(double solute) {
		return Math.pow(solute,_h)/(_KsPowH+Math.pow(solute,_h));
	}

	public double kineticDiff(double solute) {
		return _KsH*Math.pow(solute, _h-1)/(ExtraMath.sq(_KsPowH+Math.pow(solute, _h)));
	}

	public double kineticDiff(double solute, double[] paramTable, int index) {
		return Math.pow(paramTable[index], paramTable[index+1])
		        *paramTable[index+1]
		        *Math.pow(solute, paramTable[index+1]-1)
		        /(ExtraMath.sq(Math.pow(paramTable[index], paramTable[index+1])
		                +Math.pow(solute, paramTable[index+1])));
	}

	public double kineticMax() {
		return 1;
	}

	public double kineticMax(double[] paramTable, int index) {
		return 1;
	}
}
