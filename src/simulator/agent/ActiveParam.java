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

import java.util.Arrays;

import simulator.Simulator;
import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Extends SpeciesParam, adding parameters used to simulate an agents involvement in a reaction
 * 
 * Extends SpeciesParam, adding parameters used to simulate an agents involvement in a reaction
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of
 * Engineering Sciences and Applied Mathematics, Northwestern University (USA).
 *
 */
public class ActiveParam extends SpeciesParam 
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 *  Density of different cell compartments ; by convention the first one is the cell, the last one is the extracellular eps
	 */
	public Double[]   particleDensity;

	/**
	 * Solute yield from the reaction.
	 */
	public Double[][] soluteYield;
	
	/**
	 * Particle yield from the reaction.
	 */
	public Double[][] particleYield;
	
	/**
	 * Kinetic information for this reaction.
	 */
	public Double[][] reactionKinetic;

	/**
	 * For each defined solute the contribution of a cell when it dies
	 */
	public Double[]   lysisYield;

	/* ______________ UNIVERSAL MUTATOR _________________________________ */
	/**
	 * \brief Stores the reaction parameters for an active species and calculates the particle density of compounds involved in the reaction
	 * 
	 * Stores the reaction parameters for an active species and calculates the particle density of compounds involved in the reaction
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol file
	 */
	@Override
	public void init(Simulator aSim,
							XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		// Initialise simple parameter
		super.init(aSim, aSpeciesRoot,speciesDefaults);

		// Initialise particulateDensity table
		int nParticle = aSim.particleDic.size();
		int nReaction = aSim.reactionList.length;
		int nSolute = aSim.soluteList.length;

		particleDensity = ExtraMath.newDoubleArray(nParticle);
		Arrays.fill(particleDensity, Double.POSITIVE_INFINITY);

		// Do not initialise reactionKinetic using ExtraMath.newDoubleArray()
		// as the number of j-elements in each i-array varies. Each i-array is
		// cloned from the reaction mark up, via ActiveAgent.
		reactionKinetic = new Double[nReaction][];
		soluteYield = ExtraMath.newDoubleArray(nReaction, nSolute);		
		particleYield = ExtraMath.newDoubleArray(nReaction, nParticle);			
		
		int particleIndex;
		Double density;
		String name;
		for (XMLParser aParticle : aSpeciesRoot.getChildrenParsers("particle"))
		{
			name = aParticle.getName();
			particleIndex = aSim.getParticleIndex(name);
			density = aParticle.getParamConcn("density");
			if ( density.isNaN() )
			{
				aParticle = new 
					XMLParser(aSpeciesRoot.getElement().getParentElement());
				aParticle = new XMLParser(aParticle.
							getChildSuchAttribute("particle", "name", name));
				density = aParticle.getParamConcn("density");
			}
			particleDensity[particleIndex] = density;
		}
	}

}
