/**
 * \package simulator
 * \brief Package of classes that create a simulator object and capture simulation time.
 * 
 * Package of classes that create a simulator object and capture simulation time. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator;

import java.util.*;
import utils.LogFile;
import utils.ExtraMath;
import utils.XMLParser;
import simulator.geometry.*;

/**
 * \brief Class containing methods to create and query the simulation world (bulk and computational domains)
 * 
 * The simulation world is a description of all bulks and computation domains defined in the simulation. This class is used to create 
 * the world specified in the protocol files. The world mark-up in the protocol file collects the description of all bulks and 
 * computational domains defined in the simulation. The bulk mark-up defines a bulk solute compartment that is a source or sink for 
 * solutes involved in biofilm growth. The computationDomain mark-up defines the spatial region the biofilm will grow in. Only one 
 * world may be defined, but this world may contain several bulk compartments and computationDomain domains, each with a different 
 * name. Though when simulating a chemostat scenario, the name of the bulk MUST be chemostat, regardless of the corresponding 
 * computationalDomain name. 
 * 
 * This class is a component class of iDynoMiCS, released under the CECIL license. Please see www.idynomics.bham.ac.uk for more information
 * 
 * @since June 2006
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham, UK
 */
public class World 
{

	/**
	 * ArrayList of computation domains specified in a simulation scenario
	 */
	public ArrayList<Domain> domainList = new ArrayList<Domain>();
	
	/**
	 * LinkedList of bulk domains specified in a simulation scenario
	 */
	public LinkedList<Bulk>  bulkList   = new LinkedList<Bulk>();
	
	/**
	 * Holds time constraints that are applicable to this bulk
	 */
	public double[]          bulkTime   = new double[2];

	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/

	/**
	 * \brief Creates and registers the bulks and computationDomains defined in the world mark-up of the protocol file
	 * 
	 * This method, called by the createWorld method of the Simulator class, creates and registers the bulk and computation domains 
	 * that are defined in the protocol file. The computation domains are stored in an arraylist, with the bulk information stored in 
	 * a linked list
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param worldRoot	The World mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser worldRoot) 
	{
		
		try 
		{
			// Create & register the defined bulks
			for (XMLParser aBulkRoot : worldRoot.buildSetParser("bulk"))
				bulkList.add(new Bulk(aSim, aBulkRoot));
		} 
		catch(Exception e) 
		{
			LogFile.writeLog("Error trying to create bulks in World.init(): "+e);
		} 
		
		try 
		{
			// Create & register the defined domains
			for (XMLParser aDomainRoot : worldRoot.buildSetParser("computationDomain"))
				domainList.add(new Domain(aSim, aDomainRoot));
		} 
		catch(Exception e)
		{
			LogFile.writeLog("Error trying to create domains in World.init(): "+e);
		} 

	}

	/**
	 * \brief Returns a domain from the domain list that matches a specified name.
	 * 
	 * Returns a domain from the domain list that matches a specified name. One example - used to link solute grids to a particular 
	 * domain
	 *  
	 * @param cDName	The name of the domain to return from domainList	
	 * @return	The domain of the name specified by cDName, from domainList
	 */
	public Domain getDomain(String cDName) 
	{
		for (int i = 0; i<domainList.size(); i++) 
		{
			if (domainList.get(i).getName().equals(cDName)) 
			{
				return domainList.get(i);
			}
		}
		LogFile.writeLog("World.getDomain() found no domain"); 
		return null;
	}

	/**
	 * \brief Return a bulk object of the given string name
	 * 
	 * Return a bulk object of the given string name
	 * 
	 * @param bulkName	The name of the bulk that is required
	 * @return	Bulk object of that given name
	 */
	public Bulk getBulk(String bulkName) 
	{
		for (Bulk aBulk : bulkList) {
			if (aBulk.getName().equals(bulkName)) { return aBulk; }
		}
		return bulkList.getFirst();
	}

	/**
	 * \brief Determines if a bulk object of a given string name exists
	 * 
	 * Determines if a bulk object of a given string name exists
	 * 
	 * @param bulkName	The name of the bulk that is required
	 * @return	Boolean stating whether or not this bulk exists
	 */
	public boolean containsBulk(String bulkName) {
		for (Bulk aBulk : bulkList) {
			if (aBulk.getName().equals(bulkName)) { return true; }
		}
		return false;
	}

	/**
	 * \brief Calculates and returns the simulation time required to change 100% of a bulk concentration
	 * 
	 * Calculates and returns the simulation time required to change 100% of a bulk concentration
	 * 
	 * @return the time needed to change of 100% a bulk concentration
	 */
	public double getBulkTimeConstraint() {
		for(int i=1;i<bulkTime.length;i++)
			bulkTime[i]=bulkTime[i-1];

		bulkTime[0] = bulkList.getFirst().getTimeConstraint();
		for (Bulk aBulk : bulkList)
			bulkTime[0] = Math.min(bulkTime[0], aBulk.getTimeConstraint());

		return bulkTime[0];

	}

	/**
	 * \brief Returns the Sbulk value for a given solute index from that specified in the protocol file, for all bulks specified.
	 * 
	 * Returns the Sbulk value for a given solute index from that specified in the protocol file, for all bulks specified.
	 * 
	 * @param soluteIndex	The integer reference of the solute in the solute dictionary
	 * @return	An array of Sbulk values for this given solute - one for each bulk specified
	 */
	public Double[] getAllBulkValue(int soluteIndex) 
	{
		// Initialise value as a Double[] of zero's
		Double[] value = ExtraMath.newDoubleArray(bulkList.size());
		
		for (int i = 0; i < bulkList.size(); i++) 
			if (bulkList.get(i).contains(soluteIndex))
				value[i] = bulkList.get(i).getValue(soluteIndex);
		
		return value;
	}

	/**
	 * \brief Return the max value of a particular solute in the bulk
	 * 
	 * Return the max value of a particular solute in the bulk
	 * 
	 * @param soluteIndex The simulation dictionary integer reference for the solute being queried
	 * @return	Double value representing the max level of concentration in the bulk
	 */
	public double getMaxBulkValue(int soluteIndex) {
		// Rob 4/3/11: simplified
		return ExtraMath.max(getAllBulkValue(soluteIndex));
		
	}

}
