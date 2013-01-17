/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 */

/**
 * ______________________________________________________
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.agent.zoo;
import simulator.Simulator;
import simulator.agent.LocatedParam;
import utils.XMLParser;

/** Parameters common to all instances of a same species */
public class ParticulateEPSParam extends LocatedParam {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	public double transferRadius;

	public ParticulateEPSParam() {
		super();
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		double value;
		super.init(aSim, aSpeciesRoot);

		value = aSpeciesRoot.getParamLength("transferRadius");
		if(!Double.isNaN(value)) transferRadius = value;

	}
}