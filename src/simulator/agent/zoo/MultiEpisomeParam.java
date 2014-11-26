package simulator.agent.zoo;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.reaction.Reaction;
import utils.XMLParser;

public class MultiEpisomeParam extends ActiveParam
{
	LinkedList<Reaction>      pathwayKnown        = new LinkedList<Reaction>();

	public Boolean	isHighCopy = false;
	public int		nCopy, nCopyMin, nCopyMax;
	public Double	pilusLength = 0.0;
	
	public Double	exchangeLag = 0.0;
	public Double	receptionLag = 0.0;

	
	public Double	replicationSpeed = 0.26;
	
	public Double	lossProbability = 0.0;
	public Double	hotLag = 2.0;
	
	//sonia: info used in output writing
	public String	plasmidName;

	// plasmid markers (used in host compatibility)
	public ArrayList<String> hostCompatibilityMarkers = new ArrayList<String>();
	
	//plasmid Markers for plasmid compatibility (based on incompatibility groups)
	public ArrayList<String> plasmidCompatibilityMarkers = new ArrayList<String>();
	
	//sonia: transferProb is the parameter related to the plasmid transfer rate observed
	// for this plasmid.
	
	public Double transferProb = 0.0;
	//sonia 11.10.2010 - scanspeed should be a plasmid parameter 
	public Double	scanSpeed = 0.0;
	
	//sonia: 21-05-09
	//fitness cost management
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

	public HashMap<String, double[]> getSoluteYield()
	{
		// TODO Auto-generated method stub
		return null;
	}

	/* ______________ UNIVERSAL MUTATOR _________________________________ */
	/**
	 * Rob Nov 2014: I have no idea what's going on here
	 * 
	 * @param paramName
	 * @param value
	 */
	public void callParamMutator(String paramName, Double value)
	{
		try
		{
			Class[] paramTypes = null;
			paramTypes[0] = (new Double(value)).getClass();
			Method m = this.getClass().getMethod("set"+paramName, paramTypes);
			m.invoke(this, value);
		}
		catch (Exception e)
		{
			
		}
	}
}