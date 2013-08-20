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
import utils.XMLParser;

/**
 * \brief Extends Bacterium Param for use by BactEPS species, a bacterium with a permanent hydrolytic activity of its EPS capsule
 * 
 * Extends Bacterium Param for use by BactEPS species, a bacterium with a permanent hydrolytic activity of its EPS capsule.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BactEPSParam extends BacteriumParam 
{
	/**
	 * Hydrolysis speed of capsule EPS (h-1)
	 */
	public double             kHyd           = 0.007;

	/**
	 * \brief Creates a parameter storage object for the Bacterium EPS species type
	 * 
	 * Creates a parameter storage object for the Bacterium EPS species type
	 */
	public BactEPSParam() {
		super();
	}

	/**
	 * \brief Initialises Bacterium EPS Species parameters, calling the relevant superclasses to initialise common parameters and reading in the hydrolysis speed
	 * 
	 * Initialises Bacterium EPS Species parameters, calling the relevant superclasses to initialise common parameters and reading in the hydrolysis speed
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		double value;

		value = aSpeciesRoot.getParamDbl("kHyd");
		if(!Double.isNaN(value)) kHyd = value;		
	}

}
