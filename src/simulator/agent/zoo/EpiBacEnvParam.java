package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.HashMap;

import org.jdom.Element;

import simulator.Simulator;
import utils.XMLParser;

/**
 * This class extends the EpiBacParam and adds new parameters related to
 * environmental sensitivity. 
 * 
 * @author SoniaMartins
 *
 */
public class EpiBacEnvParam extends MultiEpiBacParam
{
	/**
	 * Sensitivity markers (used in environment consequences).
	 */
	public ArrayList<String> envSensitivity = new ArrayList<String>();
	
	/**
	 * List with the probabilities of dying according to the type of
	 * environment
	 */
	public HashMap<String, Double> envProbDie = new HashMap<String, Double>();
	
	/**
	 * Constructor.
	 */
	public EpiBacEnvParam()
	{
		super();
	}
	
	/**
	 * 
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot, 
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		/*
		 * Retrieving list of environments to which this species is sensitive
		 * to and the correspondent probability of dying if under the
		 * influence of that environment.
		 */
		String s;
		Double p;
		for (Element aSpeciesMarkUp :
						aSpeciesRoot.getChildrenElements("envSensitivity"))
		{
			s = aSpeciesMarkUp.getAttributeValue("name");
			envSensitivity.add(s);
			p = aSpeciesRoot.getDblAttrOfChildSuchAttribute(
									"envSensitivity","name", s, "probDie");
			envProbDie.put(s, p);
		}
	}
}
