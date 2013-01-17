
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * DiffusionSolver is an abstract class used as parent for all diffusion_solvers 
 * you could define
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

package simulator.diffusionSolver;

import simulator.diffusionSolver.multigrid.SinglegridPressure;
import simulator.geometry.IsComputationDomain;
import simulator.Simulator;
import simulator.SoluteGrid;

import utils.XMLParser;

public class Solver_pressure extends DiffusionSolver {

	// Serial version used for the serialisation of the class
	private static final long       serialVersionUID = 1L;

	protected SoluteGrid            _bLayer, _biovolume;
	protected SinglegridPressure[]  _solute;

	protected SoluteGrid[]          allSolute, allReac, allDiffReac;
	protected Simulator             mySim;
	protected IsComputationDomain   _domain;

	
	public void init(Simulator aSimulator, XMLParser xmlRoot) {
		super.init(aSimulator, xmlRoot);

		// Create the table of solute grids
		mySim = aSimulator;
		_domain = aSimulator.world.getDomain(xmlRoot.getAttribute("domain"));

		_solute = new SinglegridPressure[1];

		_bLayer = new SoluteGrid(_domain.getBiomass());
		_bLayer.gridName = "boundaryLayer";
		
		_biovolume = new SoluteGrid(_domain.getBiomass());
		_biovolume.gridName = "deltaVolume";

		_solute[0] = new SinglegridPressure(aSimulator.getSolute("pressure"), _bLayer, 0);
		_solute[0].soluteName = "pressure";
		_solute[0]._reac = _biovolume;

	}

	public void initializeConcentrationFields() {
		// Refresh then insert here the boundary layer and the diffusivity grid
		_domain.refreshBioFilmGrids();

		// We use biomass grid as boundary layer grid
		_bLayer.setGrid(_domain.getBiomass().grid);

		// Set volume change map
		updateBioVolume();

		// Set pressure to zero
		_solute[0].setSoluteGridToBulk();

	}

	@Override
	public void solveDiffusionReaction() {
		for (int j = 0; j<50; j++) {
			_solute[0].relax();
		}

	}

	public void updateBioVolume() {
		mySim.agentGrid.fitAgentVolumeRateOnGrid(_biovolume);
	}

	public SoluteGrid getPressureGrid() {
		return _solute[0]._conc;
	}

	public SoluteGrid getVolumeRateGrid() {
		return _biovolume;
	}
}
