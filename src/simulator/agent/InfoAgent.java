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
 * \brief Abstract class used for declaring Episome and MultiEpisome objects
 * 
 * Abstract class used for declaring Episome and MultiEpisome objects
 *
 */
public abstract class InfoAgent extends SpecialisedAgent implements HasSpecies {

	/**
	 * \brief Called at each time step of the simulation to compute agent characteristics
	 * 
	 * Called at each time step of the simulation to compute agent characteristics
	 */
	public abstract void internalStep();

}
