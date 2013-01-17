/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

package simulator.agent;

import simulator.Simulator;
import utils.XMLParser;

public class SpeciesParam{


	// Used to vary the initial masses of the agents
	public double initialMassCV = .1;

	/* ___________________ CONSTRUCTOR ____________________________ */
	public SpeciesParam() {
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		double value;
		value = aSpeciesRoot.getParamDbl("initialMassCV");
		if(!Double.isNaN(value)) initialMassCV = value;
	}


}
