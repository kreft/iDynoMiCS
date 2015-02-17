/**
 * \package diffusionSolver
 * \brief Package of classes used to capture the diffusion solvers that can be defined in the protocol file
 * 
 * Package of classes used to capture the diffusion solvers that can be defined in the protocol file. Solvers are used to compute 
 * the steady-state solute profile within the computational domains.
 */
package simulator.diffusionSolver;

import simulator.diffusionSolver.multigrid.SinglegridPressure;
import simulator.geometry.IsComputationDomain;
import simulator.Simulator;
import simulator.SoluteGrid;

import utils.XMLParser;

/**
 * \brief Initialises and calculates pressure effect on the simulation domain.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz
 * Centre for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 */
public class Solver_pressure extends DiffusionSolver 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long       serialVersionUID = 1L;

	/**
	 * Solute grid storing the boundary layer
	 */
	protected SoluteGrid            _bLayer;
	
	/**
	 * Solute grid storing the biomass 
	 */
	protected SoluteGrid			_biovolume;
	
	/**
	 * Pressure grid
	 */
	protected SinglegridPressure[]  _solute;
	
	/**
	 * \brief Initialise this solver, by storing the relevant solutes,
	 * boundary layer, and biomass grids required to calculate pressure.
	 * 
	 * Initialises required pressure grids
	 * 
	 * @param aSimulator	The current simulation object.
	 * @param xmlRoot	XML tags containing relevant parameters for
	 * initialising this object.
	 */
	@Override
	public void init(Simulator aSimulator, XMLParser xmlRoot)
	{
		super.init(aSimulator, xmlRoot);
		
		_solute = new SinglegridPressure[1];
		
		_bLayer = new SoluteGrid(myDomain.getBiomass());
		_bLayer.gridName = "boundaryLayer";
		
		_biovolume = new SoluteGrid(myDomain.getBiomass());
		_biovolume.gridName = "deltaVolume";
		
		_solute[0] = new SinglegridPressure(
							aSimulator.getSolute("pressure"), _bLayer, 0.0);
		_solute[0].soluteName = "pressure";
		_solute[0]._reac = _biovolume;
	}

	/**
	 * \brief Initialise concentration fields.
	 * 
	 * Refreshes domain and biomass grids and sets pressure across the grid to
	 * zero.
	 */
	@Override
	public void initializeConcentrationFields() 
	{
		// Refresh then insert here the boundary layer and the diffusivity grid.
		myDomain.refreshBioFilmGrids();
		// We use biomass grid as boundary layer grid.
		_bLayer.setGrid(myDomain.getBiomass().grid);
		// Set volume change map.
		updateBioVolume();
		// Set pressure to zero.
		_solute[0].setSoluteGridToBulk();
	}
	
	/**
	 * \brief Performs the solving algorithm on the diffusion reaction system.  
	 */
	@Override
	public void solveDiffusionReaction()
	{
		// TODO Rob 9 Feb 2015: make this parameter (50) settable in the protocol file
		for ( int j = 0; j < 50; j++ )
			_solute[0].relax();
	}

	/**
	 * \brief Update the biovolume grid after the biomass grid has been refreshed.
	 */
	public void updateBioVolume() 
	{
		mySim.agentGrid.fitAgentVolumeRateOnGrid(_biovolume);
	}

	/**
	 * \brief Return the pressure grid.
	 * 
	 * @return	Pressure grid.
	 */
	public SoluteGrid getPressureGrid()
	{
		return _solute[0]._conc;
	}

	/**
	 * \brief Return the volume rate grid.
	 * 
	 * @return	Volume rate grid.
	 */
	public SoluteGrid getVolumeRateGrid()
	{
		return _biovolume;
	}
}
