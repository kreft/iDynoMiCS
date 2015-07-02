package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.reaction.Reaction;
import utils.XMLParser;

public class MultiEpisomeParam extends ActiveParam
{
	/**
	 * Already coded in EpisomeParam
	 */
	LinkedList<Reaction>	pathwayKnown = new LinkedList<Reaction>();
	
	/**
	 * Doesn't seem to be used
	 */
	public Boolean	isHighCopy = false;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public int nCopy;
	
	/**
	 * Doesn't seem to be used
	 */
	public int nCopyMin, nCopyMax;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double	pilusLength = 0.0;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double	exchangeLag = 0.0;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double	receptionLag = 0.0;

	/**
	 * Doesn't seem to be used
	 */
	public Double	replicationSpeed = 0.26;
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double	lossProbability = 0.0;
	
	/**
	 * Doesn't seem to be used
	 */
	public Double	hotLag = 2.0;
	
	/**
	 * Info used in output writing.
	 */
	public String	plasmidName;
	
	/**
	 * Plasmid markers (used in host compatibility)
	 */
	public ArrayList<String> hostCompatibilityMarkers = new ArrayList<String>();
	
	/**
	 * Plasmid Markers for plasmid compatibility (based on incompatibility groups)
	 */
	public ArrayList<String> plasmidCompatibilityMarkers = new ArrayList<String>();
	
	/**
	 * Already coded in EpisomeParam
	 */
	public Double transferProb = 0.0;
	
	/**
	 * Scan speed should be a plasmid parameter.
	 * 
	 * TODO This is different from Brian's code, where scanSpeed belongs to
	 * the EpiBac.
	 */
	public Double	scanSpeed = 0.0;
	
	/**
	 * Fitness cost management
	 * 
	 * Replace with Brian's approach (i.e. addtional maintenance reaction)?
	 */
	public Double initialCost = 0.0;
	public Double rateDec = 0.0; 
	public Double basalCost = 0.0;

	/**
     * Called during creation of the species
     */
	public void init(Simulator aSim, XMLParser aSpeciesRoot)
	{
		isHighCopy = aSpeciesRoot.getParamBool("isHighCopy");
		//nCopy = (int) aSpeciesRoot.getParamDbl("nCopy");
		nCopy = 1;
		nCopyMin = aSpeciesRoot.getParamInt("nCopyMin");
		nCopyMax = aSpeciesRoot.getParamInt("nCopyMax");

		pilusLength = aSpeciesRoot.getParamLength("pilusLength");
		exchangeLag = aSpeciesRoot.getParamTime("exchangeLag");
		receptionLag = aSpeciesRoot.getParamTime("receptionLag");
		replicationSpeed = aSpeciesRoot.getParamDbl("replicationSpeed");
		lossProbability = aSpeciesRoot.getParamDbl("lossProbability");
		
		transferProb = aSpeciesRoot.getParamDbl("transferProb");
	    scanSpeed = aSpeciesRoot.getParamDbl("scanSpeed");;
		
		plasmidName = aSpeciesRoot.getName();	
		
		// retrieving host markers names
		for (XMLParser aChild : aSpeciesRoot.getChildrenParsers("Marker"))
			hostCompatibilityMarkers.add(aChild.getName());
		
		// retrieving plasmid markers names
		for (XMLParser aChild : aSpeciesRoot.getChildrenParsers("Compatibility"))
			plasmidCompatibilityMarkers.add(aChild.getName());
		
		for (XMLParser aChild : aSpeciesRoot.getChildrenParsers("fitnessCost"))
		{
			initialCost = aChild.getParamDbl("initialCost");
			rateDec = aChild.getParamDbl("rateOfDecay");
			basalCost = aChild.getParamDbl("basalCost");
		}
		
		System.out.println("rate of Decay is  " + rateDec);
		//System.out.println("compatibility markers being read from xml file are " + compatibilityMarkers);

	}
}