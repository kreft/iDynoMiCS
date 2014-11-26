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
		String out = aSpeciesRoot.getParam(paramName);
		return ( out == null ) ? speciesDefaults.getParam(paramName) : out;
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
		Integer out = aSpeciesRoot.getParamInt(paramName);
		return (out == XMLParser.nullInt) ?
								speciesDefaults.getParamInt(paramName) : out;
	}

	public Double getSpeciesParameterDouble(String paramName,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double out = aSpeciesRoot.getParamDbl(paramName);
		return (out == XMLParser.nullDbl) ? 
								speciesDefaults.getParamDbl(paramName) : out;
	}

	public Double getSpeciesParameterLength(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double out = aSpeciesRoot.getParamLength(paramName);
		return (out == XMLParser.nullDbl) ? 
							speciesDefaults.getParamLength(paramName) : out;
	}
	
	public Double getSpeciesParameterMass(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double out = aSpeciesRoot.getParamMass(paramName);
		return (out == XMLParser.nullDbl) ? 
								speciesDefaults.getParamMass(paramName) : out;
	}
	
	public Double getSpeciesParameterTime(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double out = aSpeciesRoot.getParamTime(paramName);
		return (out == XMLParser.nullDbl) ? 
								speciesDefaults.getParamTime(paramName) : out;
	}
	
	public Double getSpeciesParameterConcn(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Double out = aSpeciesRoot.getParamConcn(paramName);
		return (out == XMLParser.nullDbl) ? 
							speciesDefaults.getParamConcn(paramName) : out;
	}
	
	public Boolean getSpeciesParameterBool(String paramName, 
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		Boolean out = aSpeciesRoot.getParamBool(paramName);
		return ( out == XMLParser.nullBool ) ?
								speciesDefaults.getParamBool(paramName) : out;
	}
}