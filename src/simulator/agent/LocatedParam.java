/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation and their participation in relevant reactions
 * 
 * Package of utilities that create and manage agents in the simulation and their participation in relevant reactions. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;

import simulator.Simulator;
import utils.XMLParser;

/**
 * \brief Extends ActiveParam by adding location and behaviour parameters to a
 * defined agent.
 * 
 * These parameters include the radius at which division occurs, the radius at
 * which death occurs, and for self attach simulations, parameters involved in
 * simulating a cells run and tumble motion.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class LocatedParam extends ActiveParam 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	@SuppressWarnings("unused")
	private static final long serialVersionUID = 1L;
	
	/**
	 * Division radius (in µm).
	 */
	public double divRadius = 0.97;
	
	/**
	 * Degree of stochasticity used in determining cell division behaviour.
	 */
	public double divRadiusCV = 0.1;
	
	/**
	 * Fraction of the mother's total mass that is given to the baby.
	 * 
	 * Formerly called splitRatio.
	 */
	public double babyMassFrac = 0.5;
	
	/**
	 * Degree of stochasticity used in determining cell mass distribution.
	 */
	public double babyMassFracCV = 0.1;
	
	/**
	 * Minimal radius before death (in µm).
	 */
	public double deathRadius = 0.2;
	
	/**
	 * Degree of stochasticity used in determining cell death.
	 */
	public double deathRadiusCV = 0.1;
	
	/**
	 * Multiplier of the full radius to enhance distance between cells.
	 */
	public double shoveFactor = 1.15;
	
	/**
	 * Minimal distance between two cells (after shovingRadius computation).
	 */
	public double shoveLimit = 0.0;
	
	/**
	 * For simulations that model self attachment to the substratum, the agents
	 * move from the boundary layer in a random walk. This parameter captures
	 * the speed of that move.
	 * 
	 * @author Kieran Alden 170513
	 */
	public double cellRunSpeed;
	
	/**
	 * For simulations that model self attachment to the substratum, the agents
	 * move from the boundary layer in a random walk. This parameter captures
	 * the interval at which the cell will 'tumble' and change direction.
	 * 
	 * @author Kieran Alden 170513
	 */
	public double tumbleInterval;
	
	/**
	 * Some cells (e.g. some E. coli strains) will express molecules on the
	 * surface that will stick to any other surface. For versions of the
	 * simulation where self-attachment is being modelled, this parameter
	 * captures the extra cell dimension that needs to be considered in
	 * collision detection (to see if the agent sticks).
	 * 
	 * @author Kieran Alden 170513
	 */
	public double stickinessAddition;
	
	/**
	 * \brief Create a new LocatedParam parameter storage object, calling the
	 * relevant extended class constructors.
	 */
	public LocatedParam()
	{
		super();
	}

	/**
	 * \brief Assigns values to each of the location and behavioural specific
	 * parameters, reading these from the protocol file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol file.
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		Double value;
		
		// AUGUST 2013 - Change such that these can be declared as defaults, rather than for EVERY species
		// But can be overriden for each species, so need to check if in species - if not then check the defaults
		// If not in defaults, the default value hard coded into iDynoMiCS (if present) will be used
		
		value = getSpeciesParameterLength("divRadius",aSpeciesRoot,speciesDefaults);
		divRadius = value.isNaN() ? divRadius : value;
		
		value = getSpeciesParameterDouble("divRadiusCV",aSpeciesRoot,speciesDefaults);
		divRadiusCV = value.isNaN() ? divRadiusCV : value;
		
		value = getSpeciesParameterLength("deathRadius",aSpeciesRoot,speciesDefaults);
		deathRadius = value.isNaN() ? deathRadius : value;
		
		value = getSpeciesParameterDouble("deathRadiusCV",aSpeciesRoot,speciesDefaults);
		deathRadiusCV = value.isNaN() ? deathRadiusCV : value;
		
		value = getSpeciesParameterDouble("babyMassFrac",aSpeciesRoot,speciesDefaults);
		babyMassFrac = value.isNaN() ? babyMassFrac : value;
		
		value = getSpeciesParameterDouble("babyMassFracCV",aSpeciesRoot,speciesDefaults);
		babyMassFracCV = value.isNaN() ? babyMassFracCV : value;
		
		value = getSpeciesParameterLength("shoveLimit",aSpeciesRoot,speciesDefaults);
		shoveLimit = value.isNaN() ? shoveLimit : value;
		
		value = getSpeciesParameterDouble("shoveFactor",aSpeciesRoot,speciesDefaults);
		shoveFactor = value.isNaN() ? shoveFactor : value;
		
		// Attachment parameters - KA 170513
		// TODO This should be read in as a speed, not as a length!
		value = getSpeciesParameterLength("cellRunSpeed",aSpeciesRoot,speciesDefaults);
		cellRunSpeed = value.isNaN() ? cellRunSpeed : value;
		
		value = getSpeciesParameterTime("tumbleInt",aSpeciesRoot,speciesDefaults);
		tumbleInterval = value.isNaN() ? tumbleInterval : value;
		
		value = getSpeciesParameterLength("stickinessAddition",aSpeciesRoot,speciesDefaults);
		stickinessAddition = value.isNaN() ? stickinessAddition : value;
	}
	
	/**
	 * \brief Return the cell run speed.
	 * 
	 * Used in agent self-attachment scenarios.
	 * 
	 * @return	Double value stating the stored cell run speed for this species
	 * of agent.
	 */
	public double getCellRunSpeed()
	{
		return cellRunSpeed;
	}
	
	/**
	 * \brief Return the cell stickiness radius. Used in agent self-attachment
	 * scenarios.
	 * 
	 * Captures the hypothesis that for some biological bacterial cells such as
	 * E. coli, there will be receptors on the outside of the cell that adhere
	 * to the biofilm or substratum surface, and thus this needs to be added to
	 * the cell radius.
	 * 
	 * @return	Double value stating the stored stickiness radius value for
	 * agents of this species.
	 */
	public Double getStickinessRadius()
	{
		return divRadius + stickinessAddition;
	}
}
