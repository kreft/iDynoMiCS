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
 * \brief Class used to store and set parameters specific for a species. Can be extended to add additional parameters
 * 
 * This class is used to store and set parameters specific for a species. Can be extended to add additional parameters
 */
public class SpeciesParam
{
	/**
	 * Double used to vary the initial masses of the agents. Defaults to 0.1, but can be set in protocol file
	 */
	public double initialMassCV = .1;

	/**
	 * \brief Read in specific parameters for this species, changing the default if required
	 * 
	 * Read in specific parameters for this species, changing the default if required
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) 
	{
		// August 2013 - changed as this now may be specified in the species defaults and not for each species
		// So use the new method
		initialMassCV = getSpeciesParameterDouble("initialMassCV",aSpeciesRoot,speciesDefaults,initialMassCV);
	}
	
	public double getSpeciesParameterLength(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults, Double defaultValue)
	{
		if(!Double.isNaN(aSpeciesRoot.getParamLength(paramName)))
		{
			return aSpeciesRoot.getParamLength(paramName);
		}
		else if(!Double.isNaN(speciesDefaults.getParamLength(paramName)))
		{
			return speciesDefaults.getParamLength(paramName);
		}
		else
		{
			return defaultValue;
		}
	}
    public double getSpeciesParameterLength(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
    {
        return getSpeciesParameterLength(paramName, aSpeciesRoot, speciesDefaults, Double.NaN);
    }
	
	public double getSpeciesParameterDouble(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults, Double defaultValue)
	{
		if(!Double.isNaN(aSpeciesRoot.getParamDbl(paramName)))
		{
			return aSpeciesRoot.getParamDbl(paramName);
		}
		else if(!Double.isNaN(speciesDefaults.getParamDbl(paramName)))
		{
			return speciesDefaults.getParamDbl(paramName);
		}
		else
		{
			return defaultValue;
		}
	}

    public double getSpeciesParameterDouble(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
    {
        return getSpeciesParameterDouble(paramName, aSpeciesRoot, speciesDefaults, Double.NaN);
    }
	
	public String getSpeciesParameterString(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if(aSpeciesRoot.getParam(paramName)!= null)
		{
			return aSpeciesRoot.getParam(paramName);
		}
		else if(speciesDefaults.getParam(paramName)!= null)
		{
			return speciesDefaults.getParam(paramName);
		}
		else
		{
			return null;
		}
	}


}
