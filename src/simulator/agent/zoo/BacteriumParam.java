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
import java.awt.Color;


/**
 * \brief Extends LocatedParam to create an object that can store all common parameters for Bacteria-based species
 * 
 * Extends LocatedParam to create an object that can store all common parameters for Bacteria-based species
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BacteriumParam extends LocatedParam 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Maximal fraction of eps volume before excretion 
	 */
	public double epsMax = .15;
	
	/**
	 * Colour assigned to eps when visualising results in POV-Ray
	 */
	public Color epsColor = Color.lightGray;
	
	/**
	 * How particle mass is distributed: false -> normal/gaussian, true -> exponential (power of 2)
	 */
	public Boolean distMethod = false;

	/**
	 * \brief Constructor to create a BacteriumParam parameter storage object
	 * 
	 * Constructor to create a BacteriumParam parameter storage object. Simply calls the LocatedAgent super constructor
	 */
	public BacteriumParam() 
	{
		super();
	}

	/**
	 * \brief Initialises Bacterium species parameters, calling the relevant superclasses to initialise common parameters for Bacterium derived species
	 * 
	 * Initialises Bacterium species parameters, calling the relevant superclasses to initialise common parameters for Bacterium derived species
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		double value;

		value = aSpeciesRoot.getParamDbl("epsMax");
		if(!Double.isNaN(value)) epsMax = value;

		String colorName;
		colorName= aSpeciesRoot.getParam("epsColor");
		if(colorName!=null) epsColor = utils.UnitConverter.getColor(colorName);

		Boolean boolTemp;
		boolTemp = aSpeciesRoot.getParamBool("distMethod");
		if (boolTemp!=null) distMethod = boolTemp;
		
	}

}
