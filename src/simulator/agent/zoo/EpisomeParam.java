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
	/**
	 * List of reactions that this plasmid codes for.
	 */
	LinkedList<Reaction> pathwayKnown = new LinkedList<Reaction>();
	
	/**
	 * Default copy number of this plasmid in a host.
	 */
	public int nCopy = 1;
	
	/**
	 * Length of the pilus associated with this plasmid for conjugation 
	 * (in um). Note that the pilus should reach between cell surfaces, not
	 * between cell centres.
	 */
	public Double pilusLength = 2.0;
	
	/**
	 * Index used for checking whether two plasmids are compatible (i.e. they
	 * can occupy the same host).
	 */
	public int compatibilityMarker = 1;
	
	/**
	 * After exchanging a plasmid of this species, a host needs to recover
	 * before it may exchange again. This parameter gives the time needed to
	 * recover (in hours). 
	 * 
	 * TODO Consider changing name to transferLag?
	 */
	public Double exchangeLag = 5.0;
	
	/**
	 * After receiving a plasmid of this species, a host needs to recover
	 * before it may receive again. This parameter gives the time needed to
	 * recover (in hours).
	 */
	public Double receptionLag = 1.0;
	
	/**
	 * Probability that this plasmid will be lost by the host during division.
	 */
	public Double lossProbability = 0.0;
	
	/**
	 * Probability that a possible transfer will be successful.
	 * 
	 * TODO Consider changing name to transferProbability?
	 */
	public Double transferProficiency = 0.0;

	/**
	 * Called during creation of the species
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		nCopy = 1;
		
		Double temp;
		
		temp = getSpeciesParameterLength("pilusLength", aSpeciesRoot,
															speciesDefaults);
		pilusLength = Double.isFinite(temp) ? temp : pilusLength;
		
		temp = getSpeciesParameterTime("exchangeLag", aSpeciesRoot,
															speciesDefaults);
		exchangeLag = Double.isFinite(temp) ? temp : exchangeLag;
		
		temp = getSpeciesParameterTime("receptionLag", aSpeciesRoot,
															speciesDefaults);
		receptionLag = Double.isFinite(temp) ? temp : receptionLag;
		
		temp = getSpeciesParameterDouble("lossProbability", aSpeciesRoot,
															speciesDefaults);
		lossProbability = Double.isFinite(temp) ? temp : lossProbability;
		
		temp = getSpeciesParameterDouble("transferProficiency", aSpeciesRoot,
															speciesDefaults);
		transferProficiency = Double.isFinite(temp)? temp:transferProficiency;
		
		Integer val = getSpeciesParameterInteger("compatibilityMarker",
											aSpeciesRoot, speciesDefaults);
		compatibilityMarker = ( val == null ) ? compatibilityMarker : val;
	}
}