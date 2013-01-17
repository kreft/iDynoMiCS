/**
 * Project iDynoMicS
 * ______________________________________________________
 * @since June 2006
 * @copyright -> see Idynomics.java
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr)
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;

import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.reaction.Reaction;
import utils.XMLParser;

public class EpisomeParam extends ActiveParam {

	LinkedList<Reaction> pathwayKnown = new LinkedList<Reaction>();

	public int nCopy;
	public double pilusLength = 2;
	public int compatibilityMarker = 1;
	public double exchangeLag = 5;
	public double receptionLag = 1;

	public double lossProbability = 0;
	public double transferProficiency = 0;

	/**
	 * Called during creation of the species
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		nCopy = 1;

		pilusLength = aSpeciesRoot.getParamLength("pilusLength");
		exchangeLag = aSpeciesRoot.getParamTime("exchangeLag");
		receptionLag = aSpeciesRoot.getParamTime("receptionLag");

		lossProbability = aSpeciesRoot.getParamDbl("lossProbability");
		transferProficiency = aSpeciesRoot.getParamDbl("transferProficiency");
		compatibilityMarker = aSpeciesRoot.getParamInt("compatibilityMarker");
	}

}
