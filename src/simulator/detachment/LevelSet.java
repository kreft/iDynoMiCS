/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
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

public abstract class LevelSet {

	private int[]                    _gridDim;

	private double                   _res;
	private LocatedGroup[]           _shovingGrid;
	// _close is a list of agent groups on the biofilm/liquid border
	private LinkedList<LocatedGroup> _close, _alive;

	private double                   INF = Double.POSITIVE_INFINITY;
	private double                   timeStep;

	/**
	 * Generic constructor called to dynamically instantiate a child class
	 * object
	 */
	public static LevelSet staticBuilder(XMLParser root, AgentContainer anAgentGrid) {
		LevelSet out = (LevelSet) root.instanceCreator("simulator.detachment");

		out.init(anAgentGrid, root);
		return out;
	}

	public void init(AgentContainer anAgentGrid, XMLParser root) {

		_gridDim = anAgentGrid.getGridDescription();
		_res = anAgentGrid.getResolution();

		_shovingGrid = anAgentGrid.getShovingGrid();
		_close = new LinkedList<LocatedGroup>();
		_alive = new LinkedList<LocatedGroup>();
	}

	/**
	 * Identifies biofilm border
	 */
	public void refreshBorder(boolean evalErosion, Simulator aSim) {

		_close.clear();
		_alive.clear();
		timeStep = SimTimer.getCurrentTimeStep();

		// go through all elements
		for (LocatedGroup aGroup : _shovingGrid) {

			if (aGroup.isOutside) {
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
	 * @return a random grid element of the shoving grid where a attaching agent
	 * could land
	 */
	public LocatedGroup getLandingPoint() {
		return _close.get((int) ExtraMath.getUniRand()*_close.size());
	}

	/**
	 * Build list of groups belonging to the carrier
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
	 * Compute erosion time for the whole biofilm
	 * 
	 */
	public void computeLevelSet(Simulator aSim) {

		LocatedGroup trial;
		int index;

		while (_close.size()>0) {
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
	 * Test erosion time of a neighbour and consider it as a future element to compute
	 * @param index
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
	 * 
	 * @param index
	 * @return an erosion speed (micrometer by hour)
	 */
	protected abstract double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim);

	/**
	 * Get the new T value for level set as based on the values of the
	 * neighbours and marks the particle for detachment if t is less than the
	 * time step
	 * 
	 * @param i
	 * @param j
	 * @param k
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
	 * Check the validity of a solution
	 * @param s:solution value (tmp)
	 * @param f:the present f (t or INF)
	 * @param t:the present t (min(tplus, tminus)
	 * @return true if solution is not valid
	 */
	private boolean solutionValid(double s, double f, double t) {
		// check validity criteria
		// Rob 15/2/2011 Changed from solutionNotValid to remove double 
		// negative and improve readability!
		return ((f==INF) ? (s<t) : (s>t));
	}

	/**
	 * compute the maximum value of roots
	 * 
	 * @param tX
	 * @param tY
	 * @param tZ
	 * @param detachmentRate
	 * @return
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

	public LinkedList<LocatedGroup> getBorder() {
		return _close;
	}

}
