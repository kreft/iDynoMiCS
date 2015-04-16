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
 * \brief Class used to store and set parameters specific for a species.
 * 
 * Can be extended to add additional parameters.
 */
public class SpeciesParam
{
	/**
	 * Double used to vary the initial masses of the agents.
	 * Defaults to 0.1, but can be set in protocol file.
	 */
	public Double initialMassCV = 0.1;

	/**
	 * \brief Read in specific parameters for this species, changing the
	 * default if required.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol
	 * file.
	 */
	public void init(Simulator aSim,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double value = getSpeciesParameterDouble("initialMassCV",
											aSpeciesRoot, speciesDefaults);
		initialMassCV = (value == XMLParser.nullDbl) ? initialMassCV : value;
	}
	
	public String getSpeciesParameterString(String paramName,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParam(paramName) : 
					speciesDefaults.getParam(paramName);
	}
	
	/**
	 * 
	 * 
	 * @param paramName
	 * @param aSpeciesRoot
	 * @param speciesDefaults
	 * @return
	 */
	public Integer getSpeciesParameterInteger(String paramName,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamInt(paramName) : 
					speciesDefaults.getParamInt(paramName);
	}

	public Double getSpeciesParameterDouble(String paramName,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamDbl(paramName) : 
					speciesDefaults.getParamDbl(paramName);
	}

	public Double getSpeciesParameterLength(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamLength(paramName) : 
					speciesDefaults.getParamLength(paramName);
	}
	
	public Double getSpeciesParameterMass(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamMass(paramName) : 
					speciesDefaults.getParamMass(paramName);
	}
	
	public Double getSpeciesParameterTime(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamTime(paramName) : 
					speciesDefaults.getParamTime(paramName);
	}
	
	public Double getSpeciesParameterConcn(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamConcn(paramName) : 
					speciesDefaults.getParamConcn(paramName);
	}
	
	public Boolean getSpeciesParameterBool(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		return (aSpeciesRoot.isParamGiven(paramName)) ? 
					aSpeciesRoot.getParamBool(paramName) : 
					speciesDefaults.getParamBool(paramName);
	}
}