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

import utils.XMLParser;

/**
 * \brief Model a reaction using Linear Kinetics
 * 
 * Model a reaction using Linear Kinetics.
 *
 */
public class LinearKinetic extends IsKineticFactor 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * K
	 */
	private Double _K;
	
	/**
	 * \brief Initialise the kinetic, reading in kinetic parameter information
	 * from the protocol file and calculating any auxillaries needed for easing
	 * the kinetic calculation.
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file.
	 */
	public void init(Element defMarkUp)
	{
		_K = (new XMLParser(defMarkUp)).getParamDbl("K");
		nParam = 1;
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
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("K");
	}
	
	/**
	 * \brief Calculate the value of the kinetic from a given level of solute,
	 * an array containing parameters relating to the reaction, and an index to
	 * this array.
	 * 
	 * @param solute Double stating the solute concentration.
	 * @param paramTable Array of parameters relating to this reaction.
	 * @param index	An index to the parameter array.
	 * @return Double value stating the value of the kinetic for this
	 * concentration of solute.
	 */
	public Double kineticValue(Double solute, Double[] paramTable, int index)
	{
		return paramTable[index]*solute;
	}
	
	/**
	 * \brief Calculate the kinetic rate for a given solute concentration.
	 * 
	 * @param solute Double stating the solute concentration.
	 * @return Double value stating the kinetic rate for this solute
	 * concentration.
	 */
	public Double kineticValue(Double solute)
	{
		return _K*solute;
	}
	
	/**
	 * \brief Used to compute marginal difference kinetic values for a given
	 * solute concentration.
	 * 
	 * @param solute Solute concentration.
	 * @param paramTable Array of parameters relating to this reaction.
	 * @param index	An index to the parameter array.
	 * @return Differential of the reaction kinetic.
	 */
	public Double kineticDiff(Double solute, Double[] paramTable, int index)
	{
		return paramTable[index];
	}

	/**
	 * \brief Used to compute marginal difference kinetic values for a given
	 * solute concentration.
	 * 
	 * @param solute Solute concentration.
	 * @return Differential of the reaction kinetic.
	 */
	public Double kineticDiff(Double solute)
	{
		return _K;
	}

	
	
}
