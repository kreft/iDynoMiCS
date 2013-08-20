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

import org.jdom.Element;

import simulator.Simulator;
import utils.XMLParser;

/**
 * \brief Extends SpeciesParam, adding parameters used to simulate an agents involvement in a reaction
 * 
 * Extends SpeciesParam, adding parameters used to simulate an agents involvement in a reaction
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
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
	public double[]   particleDensity;

	/**
	 * Solute yield from the reaction
	 */
	public double[][] soluteYield;
	
	/**
	 * Particle yield from the reaction
	 */
	public double[][] particleYield;
	
	/**
	 * Kinetic information for this reaction
	 */
	public double[][] reactionKinetic;

	/**
	 * For each defined solute the contribution of a cell when it dies
	 */
	public double[]   lysisYield;

	/* ______________ UNIVERSAL MUTATOR _________________________________ */
	/**
	 * \brief Stores the reaction parameters for an active species and calculates the particle density of compounds involved in the reaction
	 * 
	 * Stores the reaction parameters for an active species and calculates the particle density of compounds involved in the reaction
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpeciesRoot	A Species mark-up within the specified protocol file
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) 
	{
		// Initialize simple parameter
		super.init(aSim, aSpeciesRoot,speciesDefaults);

		// Initialize particulateDensity table
		int nParticle = aSim.particleDic.size();
		int nReaction = aSim.reactionList.length;
		int nSolute = aSim.soluteList.length;

		particleDensity = new double[nParticle];
		Arrays.fill(particleDensity, Double.POSITIVE_INFINITY);

		reactionKinetic = new double[nReaction][];
		soluteYield = new double[nReaction][nSolute];		
		particleYield = new double[nReaction][nParticle];		

		XMLParser parser;
		int particleIndex;

		for (Element aChild : aSpeciesRoot.getChildren("particle")) 
		{
			// Initialize the xml parser
			parser = new XMLParser(aChild);

			// Set the density of the particular compound
			particleIndex = aSim.getParticleIndex(parser.getAttribute("name"));
			particleDensity[particleIndex] = parser.getParamConc("density");
			if(Double.isNaN(particleDensity[particleIndex])){
				// The density is specified elsewhere
				parser = new XMLParser(aSpeciesRoot.getElement().getParentElement());
				parser = new XMLParser(parser.getChildSuchAttribute("particle", "name", aSim.particleDic.get(particleIndex)));
				particleDensity[particleIndex] = parser.getParamConc("density");				
			}

		}
	}

}
