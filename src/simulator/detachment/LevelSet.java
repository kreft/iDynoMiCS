/**
 * \package simulator.detachment
 * 
 * \brief Package of classes that capture detachment of agents from the biomass
 * 
 * Package of classes that capture detachment of agents from the biomass. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.detachment;

import idyno.SimTimer;

import java.util.LinkedList;
import java.util.Collections;

import simulator.AgentContainer;
import simulator.Simulator;

import simulator.agent.LocatedGroup;
import utils.ExtraMath;
import utils.XMLParser;

/**
 * \brief Solver used for modelling detachment
 * 
 * Solver used for modelling detachment
 *  
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 *
 */
public abstract class LevelSet {

	/**
	 * Grid dimensions of the associated agent grid
	 */
	private int[]                    _gridDim;

	/**
	 * Resolution of the associated agent grid
	 */
	private double                   _res;
	
	/**
	 * Shoving grid associated with the agent grid associated with this object
	 */
	public LocatedGroup[]           _shovingGrid;
	
	/**
	 * List of agent groups on the biofilm/liquid border
	 */
	private LinkedList<LocatedGroup> _close;
	
	private LinkedList<LocatedGroup> _alive;

	private double                   INF = Double.POSITIVE_INFINITY;
	private double                   timeStep;

	/**
	 * \brief Generic constructor called to dynamically instantiate a child class object
	 * 
	 * Generic constructor called to dynamically instantiate a child class object
	 * 
	 * @param root	XML tags that contain parameters relating to detachment
	 * @param anAgentGrid	Agent grid associated with this solver
	 */
	public static LevelSet staticBuilder(XMLParser root, AgentContainer anAgentGrid) 
	{
		LevelSet out = (LevelSet) root.instanceCreator("simulator.detachment");

		out.init(anAgentGrid, root);
		return out;
	}

	/**
	 * \brief Initialise this LevelSet object by taking information from the associated grid and protocol file
	 * 
	 * Initialise this LevelSet object by taking information from the associated grid and protocol file
	 * 
	 * @param anAgentGrid	Agent grid which this solver is associated to
	 * @param root	XML tags that contain parameters relating to detachment
	 */
	public void init(AgentContainer anAgentGrid, XMLParser root) {

		_gridDim = anAgentGrid.getGridDescription();
		_res = anAgentGrid.getResolution();

		_shovingGrid = anAgentGrid.getShovingGrid();
		_close = new LinkedList<LocatedGroup>();
		_alive = new LinkedList<LocatedGroup>();
	}

	/**
	 * \brief Identify the biofilm border
	 * 
	 * Identify the biofilm border
	 * 
	 * @param evalErosion	Boolean noting whether biofilm erosion should be considered
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 */
	public void refreshBorder(boolean evalErosion, Simulator aSim) 
	{
		_close.clear();
		_alive.clear();
		timeStep = SimTimer.getCurrentTimeStep();

		// go through all elements
		for (LocatedGroup aGroup : _shovingGrid) 
		{
			if (aGroup.isOutside) 
			{
				if (aGroup.status!=1) aGroup.erosionTime = INF;
				continue;
			}

			// Test if we are in carrier with NO biomass
			// (if there is biomass it should be part of the biofilm)
			if (aGroup.status==0 && aGroup.group.size()==0) {
				aGroup.erosionTime = INF;
				continue;
			}

			// Check absence of biomass
			if (aGroup.group.isEmpty()) {
				aGroup.erosionTime = 0;
				continue;
			}

			// If you reach this line, it means you are in the biofilm.

			// Compute numbers of free elements around you
			aGroup.freeNbh();


			if (aGroup.nFreeNbh==0) {
				// you are deeply into the biofilm
				aGroup.erosionTime = INF;
			} else {
				// you are on the border
				_close.add(aGroup);

				// this sets the erosion time to infinity if
				// there is no erosion (speed is zero)
				if (evalErosion) {
					double dSpeed = getLocalDetachmentSpeed(aGroup, aSim);
					if (dSpeed != 0.0)
						aGroup.erosionTime = _res / (dSpeed*aGroup.nFreeNbh);
					else
						aGroup.erosionTime = INF;
				}

			}
		}
	}

	/**
	 * \brief Return a random LocatedGroup grid element of the shoving grid where a attaching agent could land
	 * 
	 * Return a random LocatedGroup  grid element of the shoving grid where a attaching agent could land
	 * 
	 * @return	LocatedGroup from a random element of the shoving grid
	 */
	public LocatedGroup getLandingPoint()
	{
		return _close.get(ExtraMath.getUniRandInt(_close.size()));
	}

	/**
	 * \brief Build list of groups belonging to the carrier
	 * 
	 * Build list of groups belonging to the carrier
	 * 
	 */
	public void refreshCarrier() {
		_close.clear();
		_alive.clear();

		// go through all elements
		for (LocatedGroup aGroup : _shovingGrid) {
			aGroup.resetMove();

			// aGroup.updateDeltaVolume();
			if (aGroup.isCarrier) {
				_close.add(aGroup);
			}
		}
	}

	/**
	 * \brief Compute erosion time for the whole biofilm
	 * 
	 * Compute erosion time for the whole biofilm
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 */
	public void computeLevelSet(Simulator aSim) 
	{

		LocatedGroup trial;
		int index;

		while (_close.size()>0) 
		{
			// Order the elements in _close in respect to the T value
			Collections.sort(_close, new LocatedGroup.TValueComparator());
			trial = _close.removeFirst();
			_alive.addFirst(trial);

			// get all neighbours of trial that do not belong to _close
			index = trial.gridIndex;

			// X-axis neighbour
			addToCloseAndUpdate(_shovingGrid[index].nbhGroup[0][1][1], aSim);
			addToCloseAndUpdate(_shovingGrid[index].nbhGroup[2][1][1], aSim);

			// Y-axis neighbours:
			addToCloseAndUpdate(_shovingGrid[index].nbhGroup[1][0][1], aSim);
			addToCloseAndUpdate(_shovingGrid[index].nbhGroup[1][2][1], aSim);

			// Z-axis neighbours:
			if (_gridDim[2]>1) {
				addToCloseAndUpdate(_shovingGrid[index].nbhGroup[1][1][0], aSim);
				addToCloseAndUpdate(_shovingGrid[index].nbhGroup[1][1][2], aSim);
			}
		}
	}

	/**
	 * \brief Test erosion time of a neighbour and consider it as a future element to compute
	 * 
	 * Test erosion time of a neighbour and consider it as a future element to compute
	 * 
	 * @param aGroup	Located group containing neighbours to examine
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 */
	private void addToCloseAndUpdate(LocatedGroup aGroup, Simulator aSim) {
		// If this element is not in the biofilm, get out;
		if (aGroup.group.size()==0) return;
		// if (aGroup.status!=1) { return; }

		if (!_close.contains(aGroup)&!_alive.contains(aGroup)) {
			// compute the T value
			aGroup.erosionTime = computeTValue(aGroup, aSim);

			// mark all particles for detachment if it is the case
			if (aGroup.erosionTime<timeStep) aGroup.killAll();

			// add element to _close
			_close.add(aGroup);
		}
	}

	/**
	 * \brief Return the local detachment speed
	 * 
	 * Return the local detachment speed
	 * 
	 * @param aGroup	Group of located agents on the grid
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @return An erosion speed (micrometer by hour)
	 */
	protected abstract double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim);

	/**
	 * \brief Get the new T value for level set as based on the values of the neighbours and marks the particle for detachment if it is less than the time step
	 * 
	 * Get the new T value for level set as based on the values of the neighbours and marks the particle for detachment if it is less 
	 * than the time step
	 * 
	 * @param aGroup	Group of located agents to process
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file 
	 * @return next levelset value
	 */
	private double computeTValue(LocatedGroup aGroup, Simulator aSim) {
		double nbMin, nbPlus;
		double tX, tY, tZ;

		// X (does not compute the bottom neighbour for lower points)
		nbMin = aGroup.nbhGroup[0][1][1].erosionTime;
		nbPlus = aGroup.nbhGroup[2][1][1].erosionTime;
		tX = Math.min(nbMin, nbPlus);

		// Y (implements periodicity)
		nbMin = aGroup.nbhGroup[1][0][1].erosionTime;
		nbPlus = aGroup.nbhGroup[1][2][1].erosionTime;
		tY = Math.min(nbMin, nbPlus);

		// Z (implements periodicity)
		// if 2D, z values will be the same and equal to present element value
		if (_gridDim[2]==1) {
			tZ = aGroup.erosionTime;
		} else {
			nbMin = aGroup.nbhGroup[1][1][0].erosionTime;
			nbPlus = aGroup.nbhGroup[1][1][2].erosionTime;
			tZ = Math.min(nbMin, nbPlus);
		}

		// safety check if there is no erosion: all neighbors will be INF,
		// so this one should be too
		if (tX==INF && tY==INF && tZ==INF) return INF;

		// compute the solution for all possible combinations and choose the
		// maximum value among the valid solutions
		double validSolution = 0;
		double approximateSolution = 0;
		double tmp;

		// each difference value will be tried once in any combination
		for (double fi = tX;; fi = INF) {
			for (double fj = tY;; fj = INF) {
				for (double fk = tZ;; fk = INF) {
					// skip iteration where all are INF
					if (fi==INF&&fj==INF&&fk==INF) break;

					// get the roots
					tmp = computeRoots(fi, fj, fk, getLocalDetachmentSpeed(aGroup, aSim));

					// if tmp is a number, compute maximum for approximate
					// solution, else keep approximate solution
					if (!Double.isNaN(tmp)) {
						approximateSolution = Math.max(tmp, approximateSolution);

						// check if solution is valid
						if (solutionValid(tmp, fi, tX)|solutionValid(tmp, fj, tY)
								|solutionValid(tmp, fk, tZ)) {
							// if flow reaches this point, solution is valid.
							// confront with previous maximum solution
							validSolution = Math.max(tmp, validSolution);
						}
					}
					// break
					if (fk==INF) break;
				}
				// break
				if (fj==INF) break;
			}
			// break
			if (fi==INF) break;
		}

		// validity check may return invalid in special cases due to precision
		// of float computations. In these cases, the approximate solution is
		// used
		return (validSolution==0 ? approximateSolution : validSolution);
	}

	/**
	 * \brief Check the validity of a solution
	 * 
	 * Check the validity of a solution
	 * 
	 * @param s	solution value (tmp)
	 * @param f	the present f (t or INF)
	 * @param t	the present t (min(tplus, tminus)
	 * @return true if solution is not valid
	 */
	private boolean solutionValid(double s, double f, double t) {
		// check validity criteria
		// Rob 15/2/2011 Changed from solutionNotValid to remove double 
		// negative and improve readability!
		return ((f==INF) ? (s<t) : (s>t));
	}

	/**
	 * \brief Compute the maximum value of roots
	 * 
	 * Compute the maximum value of roots
	 * 
	 * @param tX	Erosion Time
	 * @param tY	Erosion Time
	 * @param tZ	Erosion Time
	 * @param detachmentRate	Local detachment speed
	 * @return	Solution of quadratic equation
	 */
	private double computeRoots(double tX, double tY, double tZ, double detachmentRate) {
		// parameters for solving quadratic equation
		double a = (tX!=INF ? 1.0f : 0)+(tY!=INF ? 1.0f : 0)+(tZ!=INF ? 1.0f : 0);

		double b = -2.0f*((tX!=INF ? tX : 0)+(tY!=INF ? tY : 0)+(tZ!=INF ? tZ : 0));

		double c = (tX!=INF ? tX*tX : 0)+(tY!=INF ? tY*tY : 0)+(tZ!=INF ? tZ*tZ : 0)
		-ExtraMath.sq(_res/detachmentRate);

		// get the 2 solutions
		double aux = Math.sqrt(b*b-4.0f*a*c);

		// Positive solution is always the valid one
		return (-b+aux)/(2.0f*a);
	}

	/**
	 * \brief Return list of agent groups on the biofilm/liquid border
	 * 
	 * Return list of agent groups on the biofilm/liquid border
	 * 
	 * @return LinkedList of agent groups on the biofilm/liquid border
	 */
	public LinkedList<LocatedGroup> getBorder() {
		return _close;
	}

}
