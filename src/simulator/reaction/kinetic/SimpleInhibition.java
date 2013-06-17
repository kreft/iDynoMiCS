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

import utils.ExtraMath;
import utils.XMLParser;

import org.jdom.Element;

/**
 * \brief Model a reaction using a simple inhibition kinetic factor
 * 
 * Model a reaction using a simple inhibition kinetic factor. The solute of interest is specified as an attribute S and the half-maximum
 * concentration (Ki) as a parameter. Formula: mu = Ki/Ki+S
 *
 */
public class SimpleInhibition extends IsKineticFactor 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Half maximum concentration of the solute
	 */
	private double _Ki;

	/**
	 * \brief Constructor to set kinetic parameters to the supplied values
	 * 
	 * Constructor to set kinetic parameters to the supplied values
	 * 
	 * @param Ki	Concentration of solute
	 */
	public SimpleInhibition(double Ki)
	{
		_Ki = Ki;
		nParam = 1;
	}
	
	/**
	 * \brief Initialise the kinetic, reading in kinetic parameter information from the protocol file and calculating any auxillaries needed for easing the kinetic calculation
	 * 
	 * Initialise the kinetic, reading in kinetic parameter information from the protocol file and calculating any auxillaries needed 
	 * for easing the kinetic calculation
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file
	 */
	public void init(Element defMarkUp){
		_Ki = (new XMLParser(defMarkUp)).getParamDbl("Ki");		
		nParam = 1;
	}

	/**
	 * \brief Initialise the reaction from a parent of the agent
	 * 
	 * Initialise the reaction from a parent of the agent
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file
	 * @param kineticParam	Array of parameters associated with this reaction
	 * @param paramIndex	An index to the parameter array
	 */
	public void initFromAgent(Element defMarkUp,double[] kineticParam,int paramIndex){		
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("Ki");	
	}	
	
	/**
	 * \brief Calculate the value of the kinetic from a given level of solute, an array containing parameters relating to the reaction, and an index to this array
	 * 
	 * Calculate the value of the kinetic from a given level of solute, an array containing parameters relating to the reaction, and an index to this array
	 * 
	 * @param solute	Double stating the level of that solute
	 * @param paramTable	Array of parameters relating to this reaction
	 * @param index	An index to the parameter array
	 * @return Double value stating the value of the kinetic for this level of solute
	 */
	public double kineticValue(double solute,double[] paramTable,int index) {
		return paramTable[index] / (paramTable[index] + solute);
	}
	
	/**
	 * \brief Calculate the value of the kinetic for a given level of solute
	 * 
	 * Calculate the value of the kinetic for a given level of solute
	 * 
	 * @param solute	Double stating the level of that solute
	 * @return Double value stating the value of the kinetic for this level of solute
	 * 
	 */
	public double kineticValue(double solute) {
		return _Ki / (_Ki + solute);
	}
	
	
	/**
	 * \brief Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * @param solute	Solute level
	 * @param paramTable	Array of parameters relating to this reaction
	 * @param index	An index to the parameter array
	 * @return	Level of the reaction kinetic
	 */
	public  double kineticDiff(double solute,double[] paramTable,int index) 
	{
		return -_Ki / ExtraMath.sq(paramTable[index] + solute);
	}
	
	/**
	 * \brief Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * @param solute	Solute level
	 * @return	Level of the reaction kinetic
	 */
	public  double kineticDiff(double solute) {
		return -_Ki / ExtraMath.sq(_Ki + solute);
	}

}
