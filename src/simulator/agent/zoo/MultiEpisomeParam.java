package simulator.agent.zoo;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import org.jdom.Element;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.reaction.Reaction;
import utils.XMLParser;

public class MultiEpisomeParam extends ActiveParam {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID    = 1L;

	LinkedList<Reaction>      pathwayKnown        = new LinkedList<Reaction>();

	public boolean            isHighCopy = false;
	public int                nCopy, nCopyMin, nCopyMax;
	public double             pilusLength         = 0;
	
	public double             exchangeLag         = 0;
	public double             receptionLag        = 0;

	
	public double             replicationSpeed    = .26;
	
	public double             lossProbability     = 0;
	public double             hotLag              = 2;
	
	//sonia: info used in output writing
	public String             plasmidName;

	// plasmid markers (used in host compatibility)
	public ArrayList<String> hostCompatibilityMarkers = new ArrayList<String>();
	
	//plasmid Markers for plasmid compatibility (based on incompatibility groups)
	public ArrayList<String> plasmidCompatibilityMarkers = new ArrayList<String>();
	
	//sonia: transferProb is the parameter related to the plasmid transfer rate observed
	// for this plasmid.
	
	public double transferProb = 0;
	//sonia 11.10.2010 - scanspeed should be a plasmid parameter 
	public double             scanSpeed=0;
	
	//sonia: 21-05-09
	//fitness cost management
	public double initialCost=0;
	public double rateDec=0; 
	public double basalCost=0;

	/**
     * Called during creation of the species
     */
	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		XMLParser parser;
		
		//double value;
		
		isHighCopy = aSpeciesRoot.getParamBool("isHighCopy");
		//nCopy = (int) aSpeciesRoot.getParamDbl("nCopy");
		nCopy=1;
		nCopyMin = aSpeciesRoot.getParamInt("nCopyMin");
		nCopyMax = aSpeciesRoot.getParamInt("nCopyMax");

		pilusLength = aSpeciesRoot.getParamLength("pilusLength");
		exchangeLag = aSpeciesRoot.getParamTime("exchangeLag");
		receptionLag = aSpeciesRoot.getParamTime("receptionLag");
		replicationSpeed = aSpeciesRoot.getParamDbl("replicationSpeed");
		lossProbability = aSpeciesRoot.getParamDbl("lossProbability");
		
		transferProb = aSpeciesRoot.getParamDbl("transferProb");
	    scanSpeed = aSpeciesRoot.getParamDbl("scanSpeed");;
		
		plasmidName = aSpeciesRoot.getAttribute("name");	
		
		// retrieving host markers names
		for (Element aChild : aSpeciesRoot.getChildren("Marker")) {
			// Initialize the xml parser
			parser = new XMLParser(aChild);		
		hostCompatibilityMarkers.add(parser.getAttributeStr("name"));
		}
		
		// retrieving plasmid markers names
		for (Element aChild : aSpeciesRoot.getChildren("Compatibility")) {
			parser = new XMLParser(aChild);		
		plasmidCompatibilityMarkers.add(parser.getAttributeStr("name"));
		}
		
		for (Element aChild : aSpeciesRoot.getChildren("fitnessCost")){
			parser = new XMLParser(aChild);
		initialCost = parser.getParamDbl("initialCost");
		rateDec = parser.getParamDbl("rateOfDecay");
		basalCost = parser.getParamDbl("basalCost");
		}
		
		System.out.println("rate of Decay is  " + rateDec);
		//System.out.println("compatibility markers being read from xml file are " + compatibilityMarkers);

	}

	public HashMap<String, double[]> getSoluteYield() {
		// TODO Auto-generated method stub
		return null;
	}

	/* ______________ UNIVERSAL MUTATOR _________________________________ */
	public void callParamMutator(String paramName, double value) {
		try {
			Class[] paramTypes = null;
			paramTypes[0] = (new Double(value)).getClass();

			Method m = this.getClass().getMethod("set"+paramName, paramTypes);
			m.invoke(this, value);
		} catch (Exception e) {
			
		}
	}
}
