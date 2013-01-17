/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ___________________________________________________________________________
 * BactEPS : a bacterium with a permanent hydrolytic activity of its EPS capsule 
 * 
 */

/**
 * 
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

package simulator.agent.zoo;

import simulator.Simulator;
import utils.XMLParser;

public class BactEPSParam extends BacteriumParam {
	// Hydrolysis speed of capsule EPS (h-1)
	public double             kHyd           = 0.007;

	public BactEPSParam() {
		super();
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot){
		super.init(aSim,aSpeciesRoot);
		double value;

		value = aSpeciesRoot.getParamDbl("kHyd");
		if(!Double.isNaN(value)) kHyd = value;		
	}

}
