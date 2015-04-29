/**
 * \package simulator.detachment
 * 
 * \brief Package of classes that capture detachment of agents from the
 * biofilm.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.detachment;

import idyno.SimTimer;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Collections;

import simulator.AgentContainer;
import simulator.Simulator;
import simulator.agent.LocatedGroup;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Solver used for modelling detachment.
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA).
 */
public abstract class LevelSet
{
	/**
	 * Whether the agent grid is 3D (true) or 2D (false).
	 */
	private Boolean	_is3D;
	
	/**
	 * If this is set to true and the maximum threshold (MaxTh) is crossed,
	 * the simulation will be terminated at the end of the timestep.
	 */
	private Boolean _endSimWhenMaxThCrossed = false;
	
	/**
	 * Resolution of the associated agent grid.
	 */
	private Double	_res;
	
	/**
	 * Shoving grid associated with the agent grid associated with this object.
	 */
	public LocatedGroup[]	_shovingGrid;
	
	/**
	 * List of agent groups on the biofilm/liquid border.
	 */
	private LinkedList<LocatedGroup> _close;
	
	/**
	 * 
	 */
	private LinkedList<LocatedGroup> _alive;
	
	/**
	 * Constant parameter used to determine the strength of detachment.
	 */
	protected Double kDet;
	
	/**
	 * Maximum thickness that the biofilm may reach.
	 */
	protected Double maxTh;
	
	/**
	 * 
	 */
	private Double timeStep;
	
	/**
	 * \brief Generic constructor called to dynamically instantiate a child
	 * class object.
	 * 
	 * @param root	XML tags that contain parameters relating to detachment
	 * @param anAgentGrid	Agent grid associated with this solver
	 */
	public static LevelSet staticBuilder(XMLParser root, 
												AgentContainer anAgentGrid) 
	{
		LevelSet out = (LevelSet) root.instanceCreator("simulator.detachment");
		out.init(anAgentGrid, root);
		return out;
	}
	
	/**
	 * \brief Initialise this LevelSet object by taking information from the
	 * associated grid and protocol file.
	 * 
	 * @param anAgentGrid	Agent grid which this solver is associated to.
	 * @param root	XML tags that contain parameters relating to detachment.
	 */
	public void init(AgentContainer anAgentGrid, XMLParser root)
	{
		_is3D = anAgentGrid.is3D;
		_res = anAgentGrid.getResolution();
		_shovingGrid = anAgentGrid.getShovingGrid();
		_close = new LinkedList<LocatedGroup>();
		_alive = new LinkedList<LocatedGroup>();
		
		kDet = root.getParamDbl("kDet");
		Double value = root.getParamDbl("maxTh");
		maxTh = (Double.isNaN(value)? Double.POSITIVE_INFINITY:value);
		
		Boolean bool = root.getParamBool("endSimWhenMaxThCrossed");
		if ( bool != XMLParser.nullBool )
			_endSimWhenMaxThCrossed = bool;
	}
	
	/**
	 * \brief Identify the biofilm border.
	 * 
	 * @param evalErosion Boolean noting whether biofilm erosion should be
	 * considered.
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 */
	public void refreshBorder(boolean evalErosion, Simulator aSim) 
	{
		_close.clear();
		_alive.clear();
		timeStep = SimTimer.getCurrentTimeStep();
		// Go through all elements.
		for (LocatedGroup aGroup : _shovingGrid) 
		{
			if (aGroup.isOutside) 
			{
				if ( aGroup.status != 1 )
					aGroup.erosionTime = Double.POSITIVE_INFINITY;
				continue;
			}
			// Test if we are in carrier with NO biomass
			// TODO These next 2 look very similar!
			// (if there is biomass it should be part of the biofilm)
			if ( aGroup.status == 0 && aGroup.group.size() == 0 )
			{
				aGroup.erosionTime = Double.POSITIVE_INFINITY;
				continue;
			}
			// Check absence of biomass.
			if ( aGroup.group.isEmpty() )
			{
				aGroup.erosionTime = 0.0;
				continue;
			}
			// If you reach this line, it means you are in the biofilm.
			// Compute numbers of free elements around you.
			if ( aGroup.freeNbh() == 0 )
			{
				// you are deeply into the biofilm
				aGroup.erosionTime = Double.POSITIVE_INFINITY;
			}
			else
			{
				// You are on the border.
				_close.add(aGroup); 
				/*
				 * This sets the erosion time to infinity if there is no
				 * erosion (speed is zero)
				 */
				if (evalErosion)
				{
					Double dSpeed = getLocalDetachmentSpeed(aGroup, aSim);
					if ( dSpeed.equals(0.0) )
						aGroup.erosionTime = Double.POSITIVE_INFINITY;
					else
						aGroup.erosionTime = _res / (dSpeed*aGroup.nFreeNbh);
				}
			}
		}
	}
	
	/**
	 * \brief Build list of groups belonging to the carrier.
	 */
	public void refreshCarrier()
	{
		_close.clear();
		_alive.clear();
		for (LocatedGroup aGroup : _shovingGrid)
		{
			aGroup.resetMove();
			if ( aGroup.isCarrier )
				_close.add(aGroup);
		}
	}
	
	/**
	 * \brief Compute erosion time for the whole biofilm.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 */
	public void computeLevelSet(Simulator aSim) 
	{
		LocatedGroup trial;
		int i;
		while (_close.size()>0) 
		{
			// Order the elements in _close in respect to the T value.
			Collections.sort(_close, new LocatedGroup.TValueComparator());
			trial = _close.removeFirst();
			_alive.addFirst(trial);
			// Get all neighbours of trial that do not belong to _close.
			i = trial.gridIndex;
			// X-axis neighbour
			addToCloseAndUpdate(_shovingGrid[i].nbhGroup[0][1][1], aSim);
			addToCloseAndUpdate(_shovingGrid[i].nbhGroup[2][1][1], aSim);
			// Y-axis neighbours:
			addToCloseAndUpdate(_shovingGrid[i].nbhGroup[1][0][1], aSim);
			addToCloseAndUpdate(_shovingGrid[i].nbhGroup[1][2][1], aSim);
			// Z-axis neighbours:
			if ( _is3D )
			{
				addToCloseAndUpdate(_shovingGrid[i].nbhGroup[1][1][0], aSim);
				addToCloseAndUpdate(_shovingGrid[i].nbhGroup[1][1][2], aSim);
			}
		}
	}

	/**
	 * \brief Test erosion time of a neighbour and consider it as a future
	 * element to compute.
	 * 
	 * @param aGroup	Located group containing neighbours to examine.
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 */
	private void addToCloseAndUpdate(LocatedGroup aGroup, Simulator aSim)
	{
		// If this element is not in the biofilm, get out;
		if ( aGroup.group.isEmpty() )
			return;
		
		if ( ! (_close.contains(aGroup) || _alive.contains(aGroup)) )
		{
			// Compute the T value.
			aGroup.erosionTime = computeTValue(aGroup, aSim);
			// Mark all particles for detachment if it is the case.
			if ( aGroup.erosionTime < timeStep )
				aGroup.killAll("detachment");
			// Add element to _close
			_close.add(aGroup);
		}
	}
	
	/**
	 * \brief Return the local detachment speed.
	 * 
	 * @param aGroup	Group of located agents on the grid.
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @return An erosion speed (micrometer per hour).
	 */
	protected Double getLocalDetachmentSpeed(LocatedGroup aGroup,
															Simulator aSim)
	{
		if ( aGroup.cc.x > maxTh )
		{
			if ( _endSimWhenMaxThCrossed )
			{
				aSim.continueRunning = false;
				LogFile.writeLog("Maximum threshold "+maxTh);
				LogFile.writeLog("Simulation halted as threshold crossed at"+aGroup.cc);
			}
			return Double.MAX_VALUE;
		}
		return 0.0;
	}
	
	/**
	 * \brief Get the new T value for level set as based on the values of the
	 * neighbours and marks the particle for detachment if it is less than the
	 * time step.
	 * 
	 * @param aGroup	Group of located agents to process.
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file. 
	 * @return Next levelset value.
	 */
	private Double computeTValue(LocatedGroup aGroup, Simulator aSim)
	{
		Double nbMin, nbPlus;
		Double tX, tY, tZ;
		// X (does not compute the bottom neighbour for lower points)
		nbMin = aGroup.nbhGroup[0][1][1].erosionTime;
		nbPlus = aGroup.nbhGroup[2][1][1].erosionTime;
		tX = Math.min(nbMin, nbPlus);
		// Y (implements periodicity)
		nbMin = aGroup.nbhGroup[1][0][1].erosionTime;
		nbPlus = aGroup.nbhGroup[1][2][1].erosionTime;
		tY = Math.min(nbMin, nbPlus);
		// Z (implements periodicity)
		if ( _is3D )
		{
			nbMin = aGroup.nbhGroup[1][1][0].erosionTime;
			nbPlus = aGroup.nbhGroup[1][1][2].erosionTime;
			tZ = Math.min(nbMin, nbPlus);
		}
		else
			tZ = aGroup.erosionTime;
		/*
		 * Safety check if there is no erosion: all neighbors will be infinite,
		 * so this one should be too.
		 */
		if ( tX.isInfinite() && tY.isInfinite() && tZ.isInfinite() )
			return Double.POSITIVE_INFINITY;
		Double dSpeed = getLocalDetachmentSpeed(aGroup, aSim);
		if ( dSpeed == 0.0 )
			return Double.POSITIVE_INFINITY;
		/*
		 * Compute the solution for all possible combinations and choose the
		 * maximum value among the valid solutions.
		 */
		Double validSolution = 0.0;
		Double approxSolution = 0.0;
		Double tmp = 0.0;
		// Each difference value will be tried once in any combination.
		for (Double fi : Arrays.asList(tX, Double.POSITIVE_INFINITY) )
			for (Double fj : Arrays.asList(tY, Double.POSITIVE_INFINITY) )
				for (Double fk : Arrays.asList(tZ, Double.POSITIVE_INFINITY) )
				{
					// Get the root of this quadratic.
					tmp = computeRoot(fi, fj, fk, dSpeed);
					/*
					 * If tmp is a number, compute maximum for approximate
					 * solution, else keep approximate solution.
					 */
					if ( Double.isFinite(tmp) )
					{
						approxSolution = Math.max(tmp, approxSolution);
						// Check if solution is valid.
						if ( solutionValid(tmp, fi, tX) ||
								solutionValid(tmp, fj, tY) ||
								solutionValid(tmp, fk, tZ) )
						{
							/*
							 * If flow reaches this point, solution is valid.
							 * Confront with previous maximum solution.
							 */
							validSolution = Math.max(tmp, validSolution);
						}
					}
				}
		/*
		 * Validity check may return invalid in special cases due to precision
		 * of float computations. In these cases, the approximate solution is
		 * used.
		 */
		return ( validSolution.equals(0.0) ) ? approxSolution : validSolution;
	}
	
	/**
	 * \brief Check the validity of a solution.
	 * 
	 * @param s	solution value (tmp)
	 * @param f	the present f (t or INF)
	 * @param t	the present t (min(tplus, tminus)
	 * @return true if solution is not valid
	 */
	private Boolean solutionValid(Double s, Double f, Double t)
	{
		return ( f.isInfinite() ) ? ( s < t ) : ( s > t );
	}
	
	/**
	 * \brief Compute the maximum value of roots
	 * 
	 * @param tX	Erosion Time
	 * @param tY	Erosion Time
	 * @param tZ	Erosion Time
	 * @param detachmentRate	Local detachment speed
	 * @return	Solution of quadratic equation
	 */
	private Double computeRoot(Double tX, Double tY, Double tZ,
														Double detachmentRate)
	{
		// Parameters for solving quadratic equation.
		Double a = 0.0, b = 0.0, c = -ExtraMath.sq(_res/detachmentRate);;
		for (Double tVal : Arrays.asList(tX, tY, tZ) )
			if ( Double.isFinite(tVal) )
			{
				a += 1.0;
				b -= 2.0*tVal;
				c += ExtraMath.sq(tVal);
			}
		// If all tVals are infinite, return infinity
		if ( a.equals(0.0) )
			return Double.POSITIVE_INFINITY;
		// Get the 2 solutions.
		Double aux = Math.sqrt(ExtraMath.sq(b) - 4.0*a*c);
		// Positive solution is always the valid one
		return (-b+aux)/(2.0*a);
	}
	
	/**
	 * \brief Return list of agent groups on the biofilm/liquid border.
	 * 
	 * @return LinkedList of agent groups on the biofilm/liquid border.
	 */
	public LinkedList<LocatedGroup> getBorder()
	{
		return _close;
	}
}
