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

import java.io.FileWriter;
import java.io.IOException;
import simulator.AgentContainer;
import simulator.Simulator;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Extension of Agent class, adds location and parameter information
 * for an object of a particular species in the simulation.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public abstract class SpecialisedAgent extends Agent implements HasSpecies, Cloneable 
{
	/**
	 * Type of species that this agent is representing
	 */
	protected Species        _species;
	
	/**
	 * Integer index to that species within the simulation dictionary
	 */
	public int               speciesIndex;
	
	/**
	 * Boolean noting whether this agent is still active in the simulation
	 */
	public Boolean           isDead = false;
	
	/**
	 * Set of parameters associated with this specialised agent
	 */
	protected SpeciesParam   _speciesParam;
	
	/**
	 * Grid in which this agent is contained
	 */
	protected AgentContainer _agentGrid;
	
	/**
	 * Reason for agent's death. Added by Sonia Martins April 2010
	 */
	public String death;

	/**
	 * \brief Creates a SpecialisedAgent object and initialises the object in which associated parameters are stored
	 * 
	 * Creates a SpecialisedAgent object and initialises the object in which associated parameters are stored
	 */
	public SpecialisedAgent() 
	{
		// Call constructor of parent class
		super();
		_speciesParam = new SpeciesParam();
	}

	/**
	 * \brief Creates an agent of the specified species and notes the grid in which this is assigned
	 *
	 * Creates an agent of the specified species and notes the grid in which this is assigned
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A species mark-up within the specified protocol file
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) 
	{
		try 
		{
			// Create the agent object
			super.initFromProtocolFile(aSim, aSpeciesRoot);
			_agentGrid = aSim.agentGrid;
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Creating "+this.getSpecies().speciesName);
			System.exit(-1);
		}
	}

	/**
	 * \brief Clones this agent object, creating a new progeny of this agent. Ensures new clone inherits same parameters as parents
	 * 
	 * Clones this agent object, creating a new progeny of this agent. Ensures new clone inherits same parameters as parents
	 * 
	 * @throws CloneNotSupportedException	Exception should the class not implement Cloneable
	 */
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		SpecialisedAgent out = (SpecialisedAgent) super.clone();

		// Copy the references (superficial copy)
		out._species = this._species;
		out._speciesParam = this._speciesParam;
		return (Object) out;
	}

	/**
	 * \brief Create a new agent with mutated parameters based on species default values
	 * 
	 * Create a new agent with mutated parameters based on species default values
	 */
	@Override
	public abstract void createNewAgent();

	/**
	 * \brief Implemented by classes that extend this class - obtain another instance of the same species (totally independent)
	 * 
	 * Implemented by classes that extend this class - obtain another instance of the same species (totally independent)
	 */
	@Override
	public abstract SpecialisedAgent sendNewAgent() throws CloneNotSupportedException;
	
	/**
	 * \brief Registers a created agent into a respective container. Each agent must be referenced by one such container.
	 *  
	 * Registers a created agent into a respective container. Each agent must be referenced by one such container. In this case, the 
	 * species is registered into the agent grid
	 */
	@Override
	public void registerBirth() {
		_agentGrid = _species.currentSimulator.agentGrid;
		_agentGrid.registerBirth(this);
		_species.notifyBirth();
	}

	/**
	 * \brief Notifies the simulation that this agent has become too small and is then counted as dead.
	 * 
	 * Notifies the simulation that this agent has become too small and is then counted as dead. Decreases the population of this species
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has access to any resources
	 */
	public void die(Boolean isStarving)
	{
		// If you are too small, you must die !
		// Decrease the population of your species
		_species.notifyDeath();
		isDead = true;
		_agentGrid.registerDeath(this);
	}

	/**
	 * \brief Used in the calculation of delta move in Agent Container class. Will investigate further why this class returns 0
	 * 
	 * Used in the calculation of delta move in Agent Container class. Will investigate further why this class returns 0
	 * 
	 * @return	0
	 */
	public Double move()
	{
		return 0.0;
	}


	/**
	 * \brief Set the species parameters to a specified Species Parameter object
	 * 
	 * Set the species parameters to a specified Species Parameter object
	 * 
	 * @param aSpeciesParam
	 */
	public void setSpeciesParam(SpeciesParam aSpeciesParam) {
		_speciesParam = aSpeciesParam;
	}

	/**
	 * \brief Returns the object containing a set of parameters associated with a particular agent (species)
	 * 
	 * Returns the object containing a set of parameters associated with a particular agent (species)
	 */
	@Override
	public SpeciesParam getSpeciesParam() 
	{
		return _speciesParam;
	}

	/**
	 * \brief Returns the species object that is represented by this agent
	 * 
	 * Returns the species object that is represented by this agent
	 * 
	 * @return	Object of the Species class that this agent is representing
	 */
	public Species getSpecies() {
		return _species;
	}

	/**
	 * \brief Set the progenitor Specialised agent to a specified species
	 * 
	 * Set the progenitor Specialised agent to a specified species. New agents of this species are created using this information
	 * 
	 * @param aSpecies	A species object to use as the progenitor
	 */
	@Override
	public void setSpecies(Species aSpecies) 
	{
		_species = aSpecies;
		speciesIndex = aSpecies.speciesIndex;

	}

	/**
	 * \brief Return the name of the species represented by this agent.
	 * 
	 * @return	Name of the species represented.
	 */
	public String getName()
	{
		return _species.speciesName;
	}
	
	/**
	 * \brief Models a mechanical interaction between two located agents. Implemented by extending classes (LocatedAgent)
	 * 
	 * Models a mechanical interaction between two located agents. Implemented by extending classes (LocatedAgent)
	 * 
	 * @param MUTUAL	Whether movement is shared between two agents or applied only to this one
	 * @return	The move to be applied once the shoving or pull calculations have been performed
	 */
	public Double interact(boolean MUTUAL)
	{
		return 0.0;
	}
	
	/**
	 * \brief Used in POV-Ray output to display this species - writes a colour definition to the passed-in file
	 * 
	 * Used in POV-Ray output to display this species. This writes a color definition to the passed-in file. Meant for later use in 
	 * macros. Note that this routine is put here and not in Species to allow derived agents to use different colors for different 
	 * states; EpiBac is one example, with different colors for donor, recipient, and transconjugant states
	 * 
	 * @param fr	POV-Ray output file where the colour definition should be applied
	 */
	public void writePOVColorDefinition(FileWriter fr) throws IOException
	{
		fr.write("#declare "+_species.speciesName+" = color rgb < ");
		fr.write((_species.color.getRed()) / 255.0 + " , ");
		fr.write((_species.color.getGreen()) / 255.0 + " , ");
		fr.write((_species.color.getBlue()) / 255.0 + " >");
		fr.write(";\n");
	}
}
