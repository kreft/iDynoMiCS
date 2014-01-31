/**
 * \package simulator.reaction.kinetic
 * \brief Package of kinetic factors that comprise the multiplicative terms that make up the reaction kinetics in iDynoMiCS. 
 * 
 * Package of kinetic factors that comprise the multiplicative terms that make up the reaction kinetics in iDynoMiCS. This 
 * package is part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free 
 * software.  You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and 
 * INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.reaction.kinetic;

import org.jdom.Element;
import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Model a reaction using Hill Kinetics
 * 
 * Model a reaction using Hill Kinetics. Specify the solute of interest as an attribute and the half-maximum concentration (Ks) 
 * and exponent (h) as parameters. Formula: mu = Sh/Khs+Sh
 *
 */
public class HillKinetic extends IsKineticFactor 
{

	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Half-Maximum concentration of the solute of interest.
	 */
	private Double _Ks;
	
	/**
	 * The exponent.
	 */
	private Double _h;
	
	/**
	 * Auxiliary used in the calculation of the kinetic.
	 */
	private Double _KsH;
	
	/**
	 * Auxiliary used in the calculation of the kinetic.
	 */
	private Double _KsPowH;
	
	/**
	 * \brief Constructs the kinetic, setting the half-maximum concentration
	 * and the exponent.
	 * 
	 * @param Ks Half-Maximum concentration of the solute of interest.
	 * @param h	The exponent being used in the calculation.
	 */
	public HillKinetic(Double Ks, Double h) 
	{
		_Ks = Ks;
		_h = h;
		_KsH = _Ks*_h;
		_KsPowH = Math.pow(_Ks, _h);
		nParam = 2;
	}
	
	/**
	 * \brief Initialise the kinetic, reading in kinetic parameter information
	 * from the protocol file and calculating any auxillaries needed for easing
	 * the kinetic calculation.
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file.
	 */
	public void init(Element defMarkUp)
	{
		_Ks = (new XMLParser(defMarkUp)).getParamDbl("Ks");
		_h = (new XMLParser(defMarkUp)).getParamDbl("h");
		_KsH = Math.pow(_Ks, _h)*_h;
		_KsPowH = Math.pow(_Ks, _h);
		nParam = 2;
	}
	
	/**
	 * \brief Initialise the reaction from a parent of the agent.
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file.
	 * @param kineticParam	Array of parameters associated with this reaction.
	 * @param paramIndex	An index to the parameter array.
	 */
	public void initFromAgent(Element defMarkUp, Double[] kineticParam, int paramIndex)
	{
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl( "Ks");
		kineticParam[paramIndex+1] = (new XMLParser(defMarkUp)).getParamDbl("h");
	}
	
	/**
	 * \brief Calculate the value of the kinetic from a given level of solute,
	 * an array containing parameters relating to the reaction, and an index to
	 * this array.
	 * 
	 * @param solute Double of the solute concentration.
	 * @param paramTable Array of parameters relating to this reaction.
	 * @param index	An index to the parameter array.
	 * @return Double value of the kinetic rate for this solute concentration.
	 */
	public Double kineticValue(Double solute, Double[] paramTable, int index)
	{
		Double numerator, denominator;
		numerator = Math.pow(solute, paramTable[index+1]);
		denominator = Math.pow(paramTable[index], paramTable[index+1]);
		denominator += Math.pow(solute, paramTable[index+1]);
		return numerator / denominator;
	}

	/**
	 * \brief Calculate the kinetic rate for a given solute concentration.
	 * 
	 * @param solute	Double giving the solute concentration.
	 * @return Double giving kinetic rate for this solute concentration.
	 */
	public Double kineticValue(Double solute) 
	{
		return Math.pow(solute,_h)/(_KsPowH+Math.pow(solute,_h));
	}
	
	/**
	 * \brief Used to compute the kinetic differential for a given solute
	 * concentration.
	 * 
	 * @param solute	Solute concentration.
	 * @return Double value of the kinetic differential for this solute
	 * concentration.
	 */
	public Double kineticDiff(Double solute)
	{
		return _KsH * Math.pow(solute, _h-1) /
							(ExtraMath.sq(_KsPowH + Math.pow(solute, _h)));
	}
	
	/**
	 * \brief Used to compute the kinetic differential for a given solute
	 * concentration.
	 * 
	 * @param solute	Solute concentration.
	 * @param paramTable	Array of parameters relating to this reaction.
	 * @param index	An index to the parameter array.
	 * @return Double value of the kinetic differential for this solute
	 * concentration.
	 */
	public Double kineticDiff(Double solute, Double[] paramTable, int index)
	{
		Double numerator = Math.pow(paramTable[index], paramTable[index+1]);
		numerator *= paramTable[index+1];
		numerator *= Math.pow(solute, paramTable[index+1]-1);
		Double denominator = Math.pow(paramTable[index], paramTable[index+1]);
		denominator += Math.pow(solute, paramTable[index+1]); 
		return numerator / ExtraMath.sq(denominator);
	}	
}
