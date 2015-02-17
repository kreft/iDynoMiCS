/**
 * \package diffusionSolver
 * \brief Package of classes used to capture the diffusion solvers that can be
 * defined in the protocol file
 * 
 * Solvers are used to compute the solute profile within the computational
 * domains. This package is part of iDynoMiCS v1.2, governed by the CeCILL 
 * license under French law and abides by the rules of distribution of free
 * software. You can use, modify and/ or redistribute iDynoMiCS under the
 * terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the
 * following URL  "http://www.cecill.info".
 */
package simulator.diffusionSolver;

import idyno.SimTimer;
import simulator.diffusionSolver.multigrid.MultigridSolute;
import simulator.geometry.Domain;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.geometry.boundaryConditions.BoundaryAgar;
import simulator.geometry.boundaryConditions.ConnectedBoundary;
import simulator.Simulator;
import simulator.SoluteGrid;
import utils.XMLParser;

/**
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of
 * Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 */
public class Solver_multigrid extends DiffusionSolver
{
	protected MultigridSolute _bLayer;
	
	protected MultigridSolute	_diffusivity;
	protected MultigridSolute[] _solute;
	
	protected MultigridSolute[]	_biomass;

	protected SoluteGrid[]      allSolute;
	
	protected SoluteGrid[]		allReac;
	
	protected SoluteGrid[]		allDiffReac;

	protected static int        iSolute;
	
	protected static int		order;
	protected int               maxOrder;
	
	/**
	 * Number of solutes SOLVED by THIS solver
	 */
	protected int				nSolute;
	
	/**
	 * 
	 */
	protected int				nReaction;
	
	/**
	 * Number of times to relax the coarsest grid. Set in the protocol file. 
	 */
	protected int               nCoarseStep;
	
	/**
	 * Number of V-cycles to perform during each time step. Set in the
	 * protocol file.
	 */
	protected int				vCycles;
	
	/**
	 * Number of times to relax each multigrid during the downward stroke of a
	 * V-cycle. Set in the protocol file.
	 */
	protected int				nPreSteps;
	
	/**
	 * Number of times to relax each multigrid during the upward stroke of a
	 * V-cycle. Set in the protocol file.
	 */
	protected int				nPostSteps;
	
	/**
	 * 
	 */
	@Override
	public void init(Simulator aSimulator, XMLParser xmlRoot)
	{
		super.init(aSimulator, xmlRoot);
		
		nCoarseStep = xmlRoot.getParamInt("coarseStep");
		vCycles = xmlRoot.getParamInt("nCycles");
		nPreSteps = xmlRoot.getParamInt("preStep");
		nPostSteps = xmlRoot.getParamInt("postStep");
		
		// Create the table of solute grids
		nSolute = _soluteList.length;
		_solute = new MultigridSolute[nSolute];
		allSolute = new SoluteGrid[nSolute];
		allReac = new SoluteGrid[nSolute];
		allDiffReac = new SoluteGrid[nSolute];
		
		_bLayer = new MultigridSolute(_soluteList[0], "boundary layer");
		_diffusivity = 
				new MultigridSolute(_soluteList[0], "relative diffusivity");
		
		Double sBulk;
		for (int i = 0; i < nSolute; i++)
		{
			if (_soluteIndex.contains(i))
			{
				sBulk = mySim.world.getMaxBulkValue(_soluteList[i].soluteIndex);
				_solute[i] = new MultigridSolute(_soluteList[i],
												_diffusivity, _bLayer, sBulk);
			}
			else
				_solute[i] = null;
		}
		/* From this moment, nSolute is the number of solutes SOLVED by THIS
		 * solver.
		 */
		nSolute = _soluteIndex.size();
		nReaction = _reactions.size();
		maxOrder = _solute[_soluteIndex.get(0)]._conc.length;
		
		// Initialize array of reactive biomasses.
		_biomass = new MultigridSolute[nReaction];
		for (int i = 0; i < nReaction; i++)
		{
			_biomass[i] = new MultigridSolute(_soluteList[0],
											_reactions.get(i).reactionName);
			_biomass[i].resetMultigridCopies(0.0);
		}
	}
	
	@Override
	public void initializeConcentrationFields()
	{
		minimalTimeStep = 0.1*SimTimer.getCurrentTimeStep();
		
		// Refresh, then insert, the boundary layer and the diffusivity grid.
		myDomain.refreshBioFilmGrids();
		
		_bLayer.setFinest(myDomain.getBoundaryLayer());
		_bLayer.restrictToCoarsest();
		_diffusivity.setFinest(myDomain.getDiffusivity());
		_diffusivity.restrictToCoarsest();
		
		// Prepare a soluteGrid with catalyst CONCENTRATION.
		for (int i = 0; i<_biomass.length; i++)
		{
			_biomass[i].resetFinest(0.0);
			_reactions.get(i).fitAgentMassOnGrid(_biomass[i].getFinest());
			_biomass[i].restrictToCoarsest();
		}

		for (int iSolute : _soluteIndex)
			_solute[iSolute].readBulk();
	}

	/**
	 * Solve by iterative relaxation.
	 */
	@Override
	public void solveDiffusionReaction()
	{
		Double timeToSolve = SimTimer.getCurrentTimeStep();
		internalIteration = 0;
		internTimeStep = timeToSolve;
		
		/* bvm note 13.7.09:
		 * This iterative loop is only passed through once because of the
		 * value of internTimeStep used above; we leave the loop as-is though
		 * to allow future use of iterates if needed.
		 */
		while ( timeToSolve > 0 )
		{
			// Compute new equilibrium concentrations.
			stepSolveDiffusionReaction();
			
			// Update bulk concentration.
			updateBulk();
			
			// Manage iterations.
			internalIteration += 1;
			timeToSolve -= internTimeStep;
		}

		// Apply results on solute grids
		for (int iSolute : _soluteIndex)
			_solute[iSolute].applyComputation();

	}

	/**
	 * One step of the solver
	 */
	public void stepSolveDiffusionReaction()
	{
		for (int iSolute : _soluteIndex)
			_solute[iSolute].resetMultigridCopies();

		// Solve chemical concentrations on coarsest grid.
		solveCoarsest();

		// Nested iteration loop.
		for (int outer = 1; outer < maxOrder; outer++)
		{
			order = outer;
			for (int iSolute : _soluteIndex)
				_solute[iSolute].initLoop(order);

			// V-cycle loop.
			for (int v = 0; v < vCycles; v++)
			{
				// Downward stroke of V.
				while ( order > 0 )
				{
					// Pre-smoothing.
					relax(nPreSteps);
					for (int iSolute : _soluteIndex)
						_solute[iSolute].downward1(order, outer);
					
					updateReacRateAndDiffRate(order-1);
					
					for (int iSolute : _soluteIndex)
						_solute[iSolute].downward2(order, outer);

					// Reduce grid value _g for good.
					order--;
				}
				
				// Bottom of V.
				solveCoarsest();
				
				// Upward stroke of V.
				while ( order < outer )
				{
					order++;
					for (int iSolute : _soluteIndex)
						_solute[iSolute].upward(order);

					for (int iSolute : _soluteIndex)
						_solute[iSolute].truncateConcToZero(order);

					// Post-smoothing.
					relax(nPostSteps);
				}

				/* Break the V-cycles if remaining error is dominated
				 * by local truncation error (see p. 884 of Numerical Recipes)
				 */
				boolean breakVCycle = true;

				updateReacRateAndDiffRate(order);
				for (int iSolute : _soluteIndex)
					breakVCycle &= _solute[iSolute].breakVCycle(order, v);

				if (breakVCycle)
					break;
			}
		}
	}

	/**
	 * \brief Update concentration in the reactor.
	 */
	public void updateBulk()
	{
		/* Update reaction rates.
		 * This yields solute change rates in fg.L-1.hr-1
		 */
		updateReacRateAndDiffRate(maxOrder-1);
		
		// Find the connected bulks and agars and update their concentration.
		for (AllBC aBC : myDomain.getAllBoundaries())
		{
			if ( aBC instanceof ConnectedBoundary )
			{
				((ConnectedBoundary) aBC).
							updateBulk(allSolute, allReac, internTimeStep);
			}
			if ( aBC instanceof BoundaryAgar )
			{
				((BoundaryAgar) aBC).
							updateAgar(allSolute, allReac, internTimeStep);
			}
		}
		
		// Refresh the bulk concentration of the multigrids.
		for (int iSolute : _soluteIndex)
			_solute[iSolute].readBulk();
	}
	
	/**
	 * Solve the coarsest grid by relaxation Coarse grid is initialised to
	 * bulk concentration.
	 */
	public void solveCoarsest()
	{
		order = 0;
		// Reset coarsest grid to bulk concentration.
		for (int iSolute : _soluteIndex)
			_solute[iSolute].setSoluteGridToBulk(order);

		// Relax NSOLVE times.
		relax(nCoarseStep);
	}

	/**
	 * Apply nIter relaxations to the grid at the current resolution.
	 * 
	 * @param nIter
	 */
	public void relax(int nIter)
	{
		for (int j = 0; j < nIter; j++)
		{
			updateReacRateAndDiffRate(order);
			for (int iSolute : _soluteIndex)
				_solute[iSolute].relax(order);
		}
	}

	/**
	 * Call all the agents and read their uptake-rate for the current
	 * concentration.
	 * 
	 * @param resOrder
	 */
	public void updateReacRateAndDiffRate(int resOrder)
	{
		// Reset rates and derivative rates grids.
		for (int iSolute : _soluteIndex)
		{
			_solute[iSolute].resetReaction(resOrder);
			allSolute[iSolute] = _solute[iSolute]._conc[resOrder];
			allReac[iSolute] = _solute[iSolute]._reac[resOrder];
			allDiffReac[iSolute] = _solute[iSolute]._diffReac[resOrder];
		}

		// Calls the agents of the guild and sums their uptake-rate
		for (int iReac = 0; iReac<_reactions.size(); iReac++)
			_reactions.get(iReac).applyReaction(allSolute, allReac,
								allDiffReac, _biomass[iReac]._conc[resOrder]);
	}

}
