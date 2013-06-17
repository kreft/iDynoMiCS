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

import java.io.Serializable;

import org.jdom.Element;


@SuppressWarnings("serial")
/**
 * \brief Interface for all kinetic reactions that can be defined in iDynoMiCS
 * 
 * Kinetic factor mark-ups comprise the multiplicative terms that make up the entire reaction kinetic. iDynoMiCS includes a suite of 
 * kinetic factor terms. This interface defines the parameters and methods that are common for them all
 * @author kieran
 *
 */
public abstract class IsKineticFactor implements Serializable
{

	/**
	 * Number of parameters that are specified as part of this reaction
	 */
	public int nParam;
	
	/**
	 * \brief Initialise the kinetic, reading in kinetic parameter information from the protocol file and calculating any auxillaries needed for easing the kinetic calculation
	 * 
	 * Initialise the kinetic, reading in kinetic parameter information from the protocol file and calculating any auxillaries needed 
	 * for easing the kinetic calculation
	 * 
	 * @param defMarkUp	XML tags that define this kinetic in the protocol file
	 */
	public abstract void init(Element defMarkUp);
	
	/**
	 * \brief Initialise the reaction from a parent of the agent
	 * 
	 * Initialise the reaction from a parent of the agent
	 * 
	 * @param aReactionRoot	XML tags that define this kinetic in the protocol file
	 * @param kineticParam	Array of parameters associated with this reaction
	 * @param paramIndex	An index to the parameter array
	 */
	public abstract void initFromAgent(Element aReactionRoot,double[] kineticParam,int paramIndex);
	
	/**
	 * \brief Calculate the value of the kinetic for a given level of solute
	 * 
	 * Calculate the value of the kinetic for a given level of solute
	 * 
	 * @param solute	Double stating the level of that solute
	 * @return Double value stating the value of the kinetic for this level of solute
	 * 
	 */
	public abstract double kineticValue(double solute);
	
	/**
	 * \brief Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * Used to compute marginal difference kinetic values for a given solute level
	 * 
	 * @param solute	Solute level
	 * @return	Level of the reaction kinetic
	 */
	public abstract double kineticDiff(double solute);
	
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
	public abstract double kineticValue(double solute,double[] paramTable,int index);
	
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
	public abstract double kineticDiff(double solute,double[] paramTable,int index);

}
