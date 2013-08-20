/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent.zoo;
import simulator.Simulator;
import simulator.agent.LocatedParam;
import utils.XMLParser;


/**
 * \brief Stores parameters common to the Particulate EPS species - representing generic extracellular polymers and contains only the 'capsule' type.
 * 
 * Stores parameters common to the Particulate EPS species - representing generic extracellular polymers and contains only the 'capsule' type.
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class ParticulateEPSParam extends LocatedParam 
{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Radius at which an agent will transfer biomass to another agent upon the agents death
	 */
	public double transferRadius;

	/**
	 * \brief Creates a parameter storage object for the Particulate EPS species type
	 * 
	 * Creates a parameter storage object for the Particulate EPS species type
	 */
	public ParticulateEPSParam() 
	{
		super();
	}

	/**
	 * \brief Initialises Particulate EPS Species parameters, calling the relevant superclasses to initialise common parameters and setting the transfer radius
	 * 
	 * Initialises Particulate EPS Species parameters, calling the relevant superclasses to initialise common parameters and setting the transfer radius
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) 
	{
		double value;
		super.init(aSim, aSpeciesRoot, speciesDefaults);

		value = aSpeciesRoot.getParamLength("transferRadius");
		if(!Double.isNaN(value)) transferRadius = value;

	}
}