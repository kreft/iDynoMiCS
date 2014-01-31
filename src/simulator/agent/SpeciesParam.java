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
	 * \brief Double used to vary the initial masses of the agents.
	 * 
	 * Defaults to 0.1, but can be set in protocol file.
	 */
	public Double initialMassCV = 0.1;
	
	/**
	 * 
	 */
	public String domainName;
	
	/**
	 * \brief Read in specific parameters for this species, changing the
	 * default if required.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol
	 * file.
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) 
	{
		Double value;
		// August 2013 - changed as this now may be specified in the species
		// defaults and not for each species
		value = getSpeciesParameterDouble("initialMassCV", aSpeciesRoot, speciesDefaults);
		initialMassCV = value.isNaN() ? initialMassCV : value;
		
		domainName = getSpeciesParameterString("computationDomain", aSpeciesRoot, speciesDefaults);
	}
	
	public String getSpeciesParameterString(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( aSpeciesRoot.getParam(paramName) !=  null )
			return aSpeciesRoot.getParam(paramName);
		if ( speciesDefaults.getParam(paramName) != null)
			return speciesDefaults.getParam(paramName);
		return null;
	}
	
	public Integer getSpeciesParameterInteger(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamInt(paramName)) )
			return aSpeciesRoot.getParamInt(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamInt(paramName)) )
			return speciesDefaults.getParamInt(paramName);
		return 0;
	}
	
	public Double getSpeciesParameterDouble(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamDbl(paramName)) )
			return aSpeciesRoot.getParamDbl(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamDbl(paramName)) )
			return speciesDefaults.getParamDbl(paramName);
		return Double.NaN;
	}
	
	public Double getSpeciesParameterLength(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamLength(paramName)) )
			return aSpeciesRoot.getParamLength(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamLength(paramName)) )
			return speciesDefaults.getParamLength(paramName);
		return Double.NaN;
	}
	
	public Double getSpeciesParameterMass(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamMass(paramName)) )
			return aSpeciesRoot.getParamMass(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamMass(paramName)) )
			return speciesDefaults.getParamMass(paramName);
		return Double.NaN;
	}
	
	public Double getSpeciesParameterTime(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamTime(paramName)) )
			return aSpeciesRoot.getParamTime(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamTime(paramName)) )
			return speciesDefaults.getParamTime(paramName);
		return Double.NaN;
	}
	
	public Double getSpeciesParameterConc(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( ! Double.isNaN(aSpeciesRoot.getParamConc(paramName)) )
			return aSpeciesRoot.getParamConc(paramName);
		if ( ! Double.isNaN(speciesDefaults.getParamConc(paramName)) )
			return speciesDefaults.getParamConc(paramName);
		return Double.NaN;
	}
	
	public Boolean getSpeciesParameterBool(String paramName, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		if ( aSpeciesRoot.getParamBool(paramName) != null )
			return aSpeciesRoot.getParamBool(paramName);
		if ( speciesDefaults.getParamBool(paramName) != null )
			return speciesDefaults.getParamBool(paramName);
		return null;
	}
}
