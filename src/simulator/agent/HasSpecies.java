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

/**
 * \brief Interface for species classes, stating the methods that these classes should include
 * 
 * Interface for species classes, stating the methods that these classes should include
 */
public interface HasSpecies 
{	
	/**
	 * \brief Create a new progenitor with random deviation of the parameters
	 * 
	 * Create a new progenitor with random deviation of the parameters
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object cannot be cloned
	 */
	public void createNewAgent() throws CloneNotSupportedException;

	/**
	 * \brief Obtain another instance of the same species (totally independent). The returned agent is NOT registered
	 * 
	 * Obtain another instance of the same species (totally independent). The returned agent is NOT registered
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object cannot be cloned
	 */
	public HasSpecies sendNewAgent() throws CloneNotSupportedException;
	
	/**
	 * \brief Set the species to a specified Species object
	 * 
	 * Set the species parameters to a specified Species object
	 * 
	 * @param aSpecies	Species object to set this Species object to
	 */
	public void setSpecies(Species aSpecies);
	
	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	public String sendHeader();
	
	/**
	 * \brief Return the set of parameters that is associated with the object of this species
	 * 
	 * Return the set of parameters that is associated with the object of this species
	 * 
	 * @return Object of BacteriumParam that stores the parameters associated with this species
	 */
	public SpeciesParam getSpeciesParam();

	/**
	 * \brief Mutates inherited parameters and distributes particle mass - either exponentially or normally, dependent on value of distMethod
	 * 
	 *  Mutates inherited parameters and distributes particle mass - either exponentially or normally, dependent on value of distMethod
	 */
	public void mutatePop();
}
