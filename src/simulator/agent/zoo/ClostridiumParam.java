package simulator.agent.zoo;

import simulator.Simulator;
import utils.XMLParser;

/**
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class ClostridiumParam extends GeneRegBacParam
{
	/**
	 * Production of Kin (in nM s-1).
	 * TODO Convert into correct units.
	 */
	Double cK = 0.1;
	
	/**
	 * Autophosphorylation (in s-1)
	 * TODO Convert into correct units (?).
	 */
	Double alpha = 0.1;
	
	public ClostridiumParam()
	{
		super();
	}
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		Double value;
		
		value = getSpeciesParameterDouble("cK", aSpeciesRoot, speciesDefaults);
		cK = Double.isNaN(value) ? cK : value;
		
		value = getSpeciesParameterDouble("alpha", aSpeciesRoot, speciesDefaults);
		alpha = Double.isNaN(value) ? alpha : value;
		
	}
	
	
	
}
