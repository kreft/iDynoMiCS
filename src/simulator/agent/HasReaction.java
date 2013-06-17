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

import simulator.reaction.Reaction;

/**
 * \brief Interface of methods that should be implemented for agent classes that are involved in solute reactions
 * 
 * Interface of methods that should be implemented for agent classes that are involved in solute reactions
 *
 */
public interface HasReaction 
{

	/**
	 * \brief Add a previously declared reaction to the list that this agent is to participate in
	 * 
	 * Add a previously declared reaction to the list that this agent is to participate in
	 * 
	 * @param aReaction	The reaction object that this agent will participate with
	 * @param useDefaultParam	Whether or not the default reaction parameters should be used in this reaction
	 */
	public void addActiveReaction(Reaction aReaction, boolean useDefaultParam);

	/**
	 * \brief Add the details of a new reaction that a particular agent is to participate in
	 * 
	 * Add the details of a new reaction that a particular agent is to participate in
	 * 
	 * @param aReaction	The reaction object that this agent will participate with
	 * @param useDefaultParam	Whether or not the default reaction parameters should be used in this reaction
	 */
	public void addReaction(Reaction aReaction, boolean useDefaultParam);

	/**
	 * \brief Remove a reaction from those that the agent is involved in
	 * 
	 * Remove a reaction from those that the agent is involved in
	 * 
	 * @param aPathway	Reaction that this agent should no longer participate in
	 */
	public void removeReaction(Reaction aPathway);

	/**
	 * \brief Switches off a particular reaction this agent is involved in but does not delete it from reaction list
	 * 
	 * Switches off a particular reaction this agent is involved in but does not delete it from reaction list
	 * 
	 * @param aPathway	The reaction object that should be switched off
	 */
	public void switchOffreaction(Reaction aPathway);

	/**
	 * \brief Switches on a particular reaction this agent is involved in. Must have been previously declared
	 * 
	 * Switches on a particular reaction this agent is involved in. Must have been previously declared
	 * 
	 * @param aPathway	The reaction object that should be switched on
	 */
	public void switchOnReaction(Reaction aPathway);

}
