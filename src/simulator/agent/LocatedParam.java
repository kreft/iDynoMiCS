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
 * \brief Extends ActiveParam by adding location and behaviour parameters to a defined agent
 * 
 * Extends ActiveParam by adding location and behaviour parameters to a defined agent. These parameters include the radius at which 
 * division occurs, the radius at which death occurs, and for self attach simulations, parameters involved in simulating a cells run 
 * and tumble motion
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class LocatedParam extends ActiveParam 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	@SuppressWarnings("unused")
	private static final long serialVersionUID = 1L;

	/**
	 * Division radius (in µm)
	 */
	public double divRadius       = .97;
	
	/**
	 * Degree of stochasticity used in determining cell division behaviour
	 */
	public double divRadiusCV     = .1;

	/**
	 * Fraction of the mother's total mass that is given to the baby. Formerly called splitRatio
	 */
	public double babyMassFrac    = .5;
	
	/**
	 * Degree of stochasticity used in determining cell mass distribution
	 */
	public double babyMassFracCV  = .1;

	/**
	 * Minimal radius before death (in m)
	 */
	public double deathRadius     = .2;
	
	/**
	 * Degree of stochasticity used in determining cell death
	 */
	public double deathRadiusCV   = .1;

	/**
	 * Multiplier of the full radius to enhance distance between cells
	 */
	public double shoveFactor     = 1.15;

	/**
	 * Minimal distance between two cells (after shovingRadius computation)
	 */
	public double shoveLimit      = 0;
	
	/**
	 * For simulations that model self attachment to the substratum, the agents move from the boundary layer in a random walk. This 
	 * parameter captures the speed of that move (KA 170513)
	 */
	public double cellRunSpeed;
	
	/**
	 * For simulations that model self attachment to the substratum, the agents move from the boundary layer in a random walk. This 
	 * parameter captures the interval at which the cell will 'tumble' and change direction (KA 170513)
	 */
	public double tumbleInterval;
	
	/**
	 * Some cells (e.g. some e-coli species) will express molecules on the surface that will stick to any other surface. For versions of 
	 * the simulation where self-attachment is being modelled, this parameter captures the extra cell dimension that needs to be considered 
	 * in collision detection (to see if the agent sticks) (KA 170513)
	 */
	public double stickinessAddition;

	/**
	 * \brief Create a new LocatedParam parameter storage object, calling the relevant extended class constructors
	 * 
	 * Create a new LocatedParam parameter storage object, calling the relevant extended class constructirs
	 */
	public LocatedParam() {
		super();
	}

	/**
	 * \brief Assigns values to each of the location and behavioural specific parameters, reading these from the protocol file
	 * 
	 * Assigns values to each of the location and behavioural specific parameters, reading these from the protocol file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) 
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);

		//sonia 28.04.2010
		//the user can define the degree of variability in the division, split and death radius
		//by defining the "parameterCV" in the protocol file

		// AUGUST 2013 - Change such that these can be declared as defaults, rather than for EVERY species
		// But can be overriden for each species, so need to check if in species - if not then check the defaults
		// If not in defaults, the default value hard coded into iDynoMiCS (if present) will be used

		divRadius = getSpeciesParameterLength("divRadius", aSpeciesRoot, speciesDefaults, divRadius);
		divRadiusCV = getSpeciesParameterDouble("divRadiusCV", aSpeciesRoot, speciesDefaults, divRadiusCV);
		deathRadius = getSpeciesParameterLength("deathRadius", aSpeciesRoot, speciesDefaults, deathRadius);
		deathRadiusCV = getSpeciesParameterDouble("deathRadiusCV", aSpeciesRoot, speciesDefaults, deathRadiusCV);
		babyMassFrac = getSpeciesParameterDouble("babyMassFrac", aSpeciesRoot, speciesDefaults, babyMassFrac);
		babyMassFracCV = getSpeciesParameterDouble("babyMassFracCV", aSpeciesRoot, speciesDefaults, babyMassFracCV);
		shoveLimit = getSpeciesParameterLength("shoveLimit", aSpeciesRoot, speciesDefaults, shoveLimit);
		shoveFactor = getSpeciesParameterLength("shoveFactor", aSpeciesRoot, speciesDefaults, shoveFactor);

		// Attachment parameters - KA 170513
		cellRunSpeed = getSpeciesParameterLength("cellRunSpeed", aSpeciesRoot, speciesDefaults, cellRunSpeed);
		tumbleInterval = getSpeciesParameterLength("tumbleInt", aSpeciesRoot, speciesDefaults, tumbleInterval);
		stickinessAddition = getSpeciesParameterLength("stickinessAddition", aSpeciesRoot, speciesDefaults, stickinessAddition);

	}
	
	/**
	 * \brief Return the cell run speed. Used in agent self-attachment scenarios
	 * 
	 * Return the cell run speed. Used in agent self-attachment scenarios
	 * 
	 * @return	Double value stating the stored cell run speed for this species of agent
	 */
	public double getCellRunSpeed()
	{
		return cellRunSpeed;
	}
	
	/**
	 * \brief Return the cell stickiness radius. Used in agent self-attachment scenarios
	 * 
	 * Return the cell stickiness radius. Used in agent self-attachment scenarios. Captures the hypothesis that for some biological 
	 * bacterial cells such as e-coli, there will be receptors on the outside of the cell that adhere to the biofilm or substratum 
	 * surface, and thus this needs to be added to the cell radius
	 * 
	 * @return	Double value stating the stored stickiness radius value for agents of this species
	 */
	public double getStickinessRadius()
	{
		return this.stickinessAddition;
	}
	
}
