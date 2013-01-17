/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Class for containing chemical solutes, that are represented by a grid. The
 * grid is padded, 3D grid
 * Diffusivity is expressed in the local time unit
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */

package simulator;

import utils.ExtraMath;
import utils.XMLParser;

import utils.UnitConverter;

import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;

public class SoluteGrid extends SpatialGrid {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	/* ______________________ PROPERTIES OF THE SOLUTE ______________________ */

	// Identifier
	public int                  soluteIndex;

	// Diffusivity in water
	public double              diffusivity;

	// Description of diffusion, carrier and bulk domains for this solute
	private Domain _domain;

	/* _____________________ CONSTRUCTOR ___________________________________ */

	/**
	 * Constructor based on the XML file
	 * @param aSim
	 * @param xmlRoot
	 */
	public SoluteGrid(Simulator aSim, XMLParser xmlRoot) {
		double value;
		StringBuffer unit = new StringBuffer("");

		gridName = xmlRoot.getAttribute("name");
		soluteIndex = aSim.getSoluteIndex(gridName);
		_domain = aSim.world.getDomain(xmlRoot.getAttribute("domain"));

		/* Set the resolution and create the grid ________________ */
		value = xmlRoot.getParamLength("resolution");
		if (Double.isNaN(value)) {
			useDomaingrid();
		} else {
			specifyResolution(value);
		}
		initGrids();
		
		/* Set the diffusivity ____________________________________ */
		value = xmlRoot.getParamDbl("diffusivity", unit);
		value *= UnitConverter.time(unit.toString());
		value *= UnitConverter.length(unit.toString());
		value *= UnitConverter.length(unit.toString());
		diffusivity = value;

		/* Set the initial concentration __________________________ */
		value = xmlRoot.getParamDbl("concentration");
		// If no value specified, use the maximal concentration of the bulks
		if (Double.isNaN(value)) value = ExtraMath.max(aSim.world.getAllBulkValue(soluteIndex));
		setAllValueAt(value);
	}


	/* The 2 next functions are used when creating the multigrids */

	public SoluteGrid(int nI, int nJ, int nK, double res) {
		super(nI, nJ, nK, res);
	}

	public SoluteGrid(int nI, int nJ, int nK, double res,String aName, Domain aDomain) {
		super(nI, nJ, nK, res);
		gridName = aName;
		_domain = aDomain;
	}

	public SoluteGrid(int nI, int nJ, int nK, double res,SoluteGrid aSolG) {
		super(nI, nJ, nK, res);
		useExternalSoluteGrid(aSolG);
	}

	public SoluteGrid(SoluteGrid aSolG){
		gridName = aSolG.gridName;
		diffusivity = aSolG.diffusivity;
		_domain = aSolG._domain;

		_reso = aSolG.getResolution();
		_nI = aSolG.getGridSizeI();
		_nJ = aSolG.getGridSizeJ();
		_nK = aSolG.getGridSizeK();

		initGrids();

	}

	public void useExternalSoluteGrid(SoluteGrid aSolG) {
		gridName = aSolG.gridName;
		soluteIndex = aSolG.soluteIndex;
		diffusivity = aSolG.diffusivity;
		_domain = aSolG._domain;
	}


	/**
	 * Use the size and the resolution used to define the computation domain to
	 * define the solute grid
	 */
	public void useDomaingrid() {
		_reso = _domain.getGrid().getResolution();
		_nI = _domain.getGrid().getGridSizeI();
		_nJ = _domain.getGrid().getGridSizeJ();
		_nK = _domain.getGrid().getGridSizeK();
	}

	/**
	 * Give size of grid for the given resolution, based on length defined in
	 * the domain
	 * 
	 * @param reso
	 */
	public void specifyResolution(double reso) {
		_reso = reso;
		_nI = (int) Math.ceil(_domain.getGrid().getGridLength(1)/_reso);
		_nJ = (int) Math.ceil(_domain.getGrid().getGridLength(2)/_reso);
		_nK = (int) Math.ceil(_domain.getGrid().getGridLength(3)/_reso);
	}

	/* ________________________ MAIN METHODS ______________________________ */

	public void refreshBoundary() {
		for (AllBC aBC:_domain.getAllBoundaries()) {
			aBC.refreshBoundary(this);
		}
	}



	/* ________________________ GET & SET __________________________________ */

	public String getName() {
		return gridName;
	}

	public double getDiffusivity() {
		return diffusivity;
	}

	public Domain getDomain() {
		return _domain;
	}

}
