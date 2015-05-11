package simulator.agent.zoo;

import simulator.Simulator;
import utils.XMLParser;

/**
 * 
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class GeneRegBacParam extends BacteriumParam 
{
	

	/**
	 * 
	 */
	public Double rtol = 0.001;
	
	/**
	 * 
	 */
	public Double hmax = 0.000001;
	
	
	
	
	public GeneRegBacParam()
	{
		super();
	}
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		Double value;
		
		value = getSpeciesParameterDouble("rtol", aSpeciesRoot, speciesDefaults);
		rtol = Double.isNaN(value) ? rtol : value;
		
		value = getSpeciesParameterDouble("hmax", aSpeciesRoot, speciesDefaults);
		hmax = Double.isNaN(value) ? hmax : value;
		
	}
	
	
	
}