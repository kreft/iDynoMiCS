package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.HashMap;

import org.jdom.Element;

import simulator.Simulator;
import utils.XMLParser;

/**
 * This class extends the EpiBacParam and adds new parameters related to environmental sensitivity. 
 * 
 * @author SoniaMartins
 *
 */
public class EpiBacEnvParam extends MultiEpiBacParam{
	

	private static final long serialVersionUID = 1L;

	//sensitivity markers (used in environment consequences...)
	public ArrayList<String> envSensitivity = new ArrayList<String>();
	
	//list with the probabilities of diying according to the type of environment
	public HashMap<String, Double> envProbDie = new HashMap<String, Double>();
	
	
	
	/**
	 * constructor
	 */
	public EpiBacEnvParam() {
		super();
	}
	
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults) {
		super.init(aSim, aSpeciesRoot, speciesDefaults);

	
	//retrieving list of environments to which this species is sensitive to and the correspondent
	//probability of dying if under the influence of that environment.
	for (Element aSpeciesMarkUp : aSpeciesRoot.buildSetMarkUp("envSensitivity")) {
		envSensitivity.add(aSpeciesMarkUp.getAttributeValue("name"));
		String s = aSpeciesMarkUp.getAttributeValue("name");
		envProbDie.put(s, aSpeciesRoot.getDblAttrOfChildSuchAttribute("envSensitivity","name", s, "probDie"));
		
	}


	}
}
