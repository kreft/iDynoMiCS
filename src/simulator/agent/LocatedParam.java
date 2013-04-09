/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 */

/**
 * 
 * @since April 2007
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * 
 */

package simulator.agent;

import simulator.Simulator;
import utils.XMLParser;

public class LocatedParam extends ActiveParam {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	/* Cell division parameters __________________________ */
	// Division radius (in µm)
	public double divRadius       = .97;
	public double divRadiusCV     = .1;

	// Rob 21/1/11: changed splitRatio to babyMassFrac, as clearer name
	// Fraction of the mother's total mass that is given to the baby
	public double babyMassFrac    = .5;
	public double babyMassFracCV  = .1;

	// Minimal radius before death (in m)
	public double deathRadius     = .2;
	public double deathRadiusCV   = .1;

	// Multiplier of the full radius to enhance distance between cells
	public double shoveFactor     = 1.15;


	// Minimal distance between two cells (after shovingRadius computation)
	public double shoveLimit      = 0;

	public LocatedParam() {
		super();
	}

	public void init(Simulator aSim, XMLParser aSpeciesRoot) {
		super.init(aSim, aSpeciesRoot);
		double value;

		//sonia 28.04.2010
		//the user can define the degree of variability in the division, split and death radius
		//by defining the "parameterCV" in the protocol file

		value = aSpeciesRoot.getParamLength("divRadius");
		if(!Double.isNaN(value)) divRadius = value;

		value = aSpeciesRoot.getParamDbl("divRadiusCV");
		if(!Double.isNaN(value)) divRadiusCV = value;

		value = aSpeciesRoot.getParamLength("deathRadius");
		if(!Double.isNaN(value)) deathRadius = value;

		value = aSpeciesRoot.getParamDbl("deathRadiusCV");
		if(!Double.isNaN(value)) deathRadiusCV = value;	

		value = aSpeciesRoot.getParamDbl("babyMassFrac");
		if(!Double.isNaN(value)) babyMassFrac  = value;

		value = aSpeciesRoot.getParamDbl("babyMassFracCV");
		if(!Double.isNaN(value)) babyMassFracCV = value;	

		value = aSpeciesRoot.getParamLength("shoveLimit");
		if(!Double.isNaN(value)) shoveLimit = value;

		value = aSpeciesRoot.getParamLength("shoveFactor");
		if(!Double.isNaN(value)) shoveFactor = value;		
	}
	
/*	public double getMaximalRadius() {
		return divRadius*(1+2*divRadiusCV);
	}
	public double getMinimalRadius() {
		return divRadius*(1-2*divRadiusCV)/2;
	}	*/
}
