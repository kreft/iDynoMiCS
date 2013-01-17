/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */

package simulator.agent;

import java.util.Arrays;

import org.jdom.Element;

import simulator.Simulator;
import utils.XMLParser;

public class ActiveParam extends SpeciesParam {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	// Density of different cell compartments ; by convention the first one is
	// the cell the last one is the extracellular eps
	public double[]   particleDensity;

	// Parameters of the reaction
	public double[][] soluteYield;
	public double[][] particleYield;
	public double[][] reactionKinetic;

	// For each defined solute the contribution of a cell when it dies
	public double[]   lysisYield;

	/* ______________ UNIVERSAL MUTATOR _________________________________ */
	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		// Initialize simple parameter
		super.init(aSim, aSpeciesRoot);

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

		for (Element aChild : aSpeciesRoot.getChildren("particle")) {
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
