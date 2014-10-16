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

import simulator.agent.ActiveParam;
import simulator.reaction.Reaction;
import simulator.Simulator;
import utils.XMLParser;

public class EpisomeParam extends ActiveParam
{

	LinkedList<Reaction> pathwayKnown = new LinkedList<Reaction>();

	public int nCopy;
	public Double pilusLength = 2.0;
	public int compatibilityMarker = 1;
	public Double exchangeLag = 5.0;
	public Double receptionLag = 1.0;

	public Double lossProbability = 0.0;
	public Double transferProficiency = 0.0;

	/**
	 * Called during creation of the species
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot)
	{
		nCopy = 1;
		
		pilusLength = aSpeciesRoot.getParamLength("pilusLength");
		exchangeLag = aSpeciesRoot.getParamTime("exchangeLag");
		receptionLag = aSpeciesRoot.getParamTime("receptionLag");
		
		lossProbability = aSpeciesRoot.getParamDbl("lossProbability");
		transferProficiency = aSpeciesRoot.getParamDbl("transferProficiency");
		compatibilityMarker = aSpeciesRoot.getParamInt("compatibilityMarker");
	}

}