/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation
 * and their participation in relevant reactions
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.agent;

import java.util.*;

import simulator.Simulator;
import simulator.AgentContainer;
import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.SoluteGrid;
import utils.ExtraMath;
import utils.LogFile;

/**
 * \brief Object to hold a group of agents in one location on the agent grid.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class LocatedGroup
{
	/**
	 * Agent container this located group is within
	 */
	public AgentContainer agentGrid;
	
	/**
	 * Linked list to hold members of this group
	 */
	public LinkedList<LocatedAgent> group = new LinkedList<LocatedAgent>();

	/**
	 * Concentration of species in this group, thus area represented
	 */
	public Double[] speciesConcentration;
	
	/**
	 * Total volume of agents in this group.
	 */
	public Double totalVolume = 0.0;
	
	/**
	 * Total concentration of agents in this group.
	 */
	public Double totalConcentration = 0.0;
	
	/**
	 * Total mass of agents in this group.
	 */
	public Double totalMass = 0.0;
	
	/**
	 * Holder for the time it would take for this group to be eroded from the
	 * top of the biofilm. 
	 */
	public Double erosionTime = Double.NaN;
	
	/**
	 * Distance of the group at this location from the carrier.
	 */
	public Double distanceFromCarrier = 0.0;
	
	/**
	 * Holder for the ratio between this group's erosionTime and the
	 * simulation timestep.
	 */
	public Double erosionRatio = 0.0;

	/**
	 * The index of the grid at which this located group represents.
	 */
	public int gridIndex;
	
	/**
	 * Coordinates of this location as a continuous vector.
	 */
	public ContinuousVector cc;
	
	/**
	 * Coordinates of this location as a discrete vector.
	 */
	public DiscreteVector dc;
	
	/**
	 * Index of the neighbours of this group.
	 */
	public int[][][] nbhIndex = new int[3][3][3];
	
	/**
	 * Neighbouring groups around this group.
	 */
	public LocatedGroup[][][] nbhGroup = new LocatedGroup[3][3][3];

	/**
	 * Space occupation
	 * -1 outside
	 * 0  carrier
	 * 1  biofilm
	 * 2  liquid
	 * 3  bulk
	 */ 
	public int status = 2;
	
	/**
	 * Boolean stating whether this location is in the carrier.
	 */
	public boolean isCarrier = false;
	
	/**
	 * Boolean stating whether this location is outside the grid.
	 */
	public boolean isOutside;
	
	/**
	 * Number of free neighbours around this location.
	 */
	public int nFreeNbh;
	
	/**
	 * Vector to hold an amount of distance an agent is to move.
	 */
	public ContinuousVector move = new ContinuousVector();

	

	/**
	 * \brief Set the coordinates of this group and check these are inside the grid
	 * 
	 *  Set the coordinates of this group and check these are inside the grid 
	 *  
	 *  @param index	Index of the grid in which this group is located
	 *  @param anAgentGrid	The agent grid associated with this index
	 *  @param aSimulator	The current simulation object 
	 */
	public LocatedGroup(int index, AgentContainer anAgentGrid, Simulator aSimulator) 
	{
		//LogFile.writeLogDebug("Debugging LocatedGroup()");
		agentGrid = anAgentGrid;
		// Spatial location of the group
		gridIndex = index;
		
		// Coordinates if padding is removed
		dc = agentGrid.getGridPosition(gridIndex);
		cc = agentGrid.getGridLocation(gridIndex);

		// Initialise biomass statistics
		speciesConcentration = ExtraMath.newDoubleArray(aSimulator.speciesDic.size());

		// Check if the group is inside the domain
		isOutside = false;
		for (AllBC aBC : anAgentGrid.domain.getAllBoundaries())
			if ( aBC.isOutside(cc) )
			{
				isOutside = true;
				status = -1;
				// setBoundary sets the status to 3!
				aBC.setBoundary(this);
				break;
			}
		/*
		if ( isOutside )
			LogFile.writeLogDebug("\tLocatedGroup at "+dc.toString()+" is outside");
		else
			LogFile.writeLogDebug("\tLocatedGroup at "+dc.toString()+" is inside");
		*/
	}

	/**
	 * \brief Builds neighbourhood reference map.
	 */
	public void init()
	{
		if ( isOutside )
			return;
		if ( agentGrid.is3D )
			testNbh_3D(agentGrid.getShovingGrid());
		else
			testNbh_2D(agentGrid.getShovingGrid());
		distanceFromBorders();
		if ( distanceFromCarrier < agentGrid.getResolution() )
			isCarrier = true;
	}
	
	/**
	 * \brief Refresh status and concentration of the group.
	 */
	public void refreshElement()
	{
		Double volume = ExtraMath.cube(agentGrid.getResolution());
		if ( ! Simulator.isChemostat)
		{
			// Refresh group status (carrier, biofilm, free)
			if ( status > 0 )
				status = ( group.size() > 0 ) ? 1 : 2;
			if ( isCarrier )
				status = 0;
		}
		
		// Refresh biomass density
		Double value = 0.0;
		totalConcentration = 0.0;
		totalMass = 0.0;
		Arrays.fill(speciesConcentration, 0.0);
		for (LocatedAgent aLoc : group)
		{
			totalMass += aLoc.getTotalMass();
			value = aLoc.getTotalMass()/volume;
			/* We treat agents as cylinders in 2D, and so the element volume
			 * is the same cube as in 3D; no concentration correction is 
			 * necessary
			 */
			totalConcentration += value;
			speciesConcentration[aLoc.speciesIndex] += value;
		}
	}

	/**
	 * \brief Refresh the volume statistics of this group.
	 * 
	 * @return	Double noting the total volume of agents in this group.
	 */
	public Double refreshVolume()
	{
		totalVolume = 0.0;
		for (LocatedAgent aLoc : group)
			totalVolume += aLoc.getVolume(true);
		return totalVolume;
	}

	/**
	 * \brief Compute the gradient due to pressure and use it to set the
	 * advective affect.
	 * 
	 * @param pressure	Pressure grid.
	 * @param deltaT	DeltaT.
	 * @return	Norm of the movement vector under the affect of the pressure
	 * gradient.
	 */
	public Double computeMove(SoluteGrid pressure, Double deltaT) 
	{
		if ( this.isOutside )
			resetMove();
		else
		{
			move = pressure.getGradient(this.cc);
			if ( move.isValid() )
				move.times(-deltaT);
			else
				resetMove();
		}
		return move.norm();
	}
	
	public void resetMove()
	{
		move.reset();
	}

	/**
	 * \brief Scale the movement vector for the grid element and apply to each
	 * agent.
	 * 
	 * @param alpha	Scaling factor to be applied to a move.
	 */
	public void addMoveToAgents(Double alpha)
	{
		for (LocatedAgent aLoc : group)
		{
			move.times(alpha);
			aLoc.addMovement(move);
		}
	}
	
	/**
	 * \brief Safely remove located agent from agentList & agentGrid.
	 * 
	 * @param reason One-word string giving the reason for death.
	 */
	public void killAll(String reason)
	{
		for ( LocatedAgent aLoc : group )
		{
			aLoc.death = reason;
			agentGrid._agentToKill.add(aLoc);
			agentGrid.agentList.remove(aLoc);
		}
		group.clear();
		if ( ! Simulator.isChemostat )
			status = 2;
	}
	
	/**
	 * \brief Remove an agent from this LocatedGroup.
	 * 
	 * @param anAgent	LocatedAgent to remove from this group.
	 */
	public void remove(LocatedAgent anAgent)
	{
		group.remove(anAgent);
		if ( group.isEmpty() && !Simulator.isChemostat )
			status = 2;
	}
	
	/**
	 * \brief Add an agent to this LocatedGroup.
	 * 
	 * @param anAgent	LocatedAgent to add to this group.
	 */
	public void add(LocatedAgent anAgent) 
	{
		group.add(anAgent);
		status = 1;
		anAgent.setGridIndex(gridIndex);
	}
	
	/**
	 * \brief Move the X parameter by a specified amount.
	 * 
	 * @param i	The current I grid element.
	 * @return	Located group reached by that move.
	 */
	public LocatedGroup moveX(int i)
	{
		int delta = Integer.signum(i);
		LocatedGroup out = nbhGroup[delta+1][1][1];
		i -= delta;
		while ( i != 0 )
		{
			out = moveX(delta);
			i -= delta;
		}
		return out;
	}
	
	/**
	 * \brief Move the Y parameter by a specified amount.
	 * 
	 * @param j	The current grid element.
	 * @return	Located group reached by that move.
	 */
	public LocatedGroup moveY(int j)
	{
		int delta = Integer.signum(j);
		LocatedGroup out = nbhGroup[1][delta+1][1];
		j -= delta;
		while ( j != 0 )
		{
			out = moveY(delta);
			j -= delta;
		}
		return out;
	}

	/**
	 * \brief Move the Z parameter by a specified amount.
	 * 
	 * @param k	The current grid element.
	 * @return	Located group reached by that move.
	 */
	public LocatedGroup moveZ(int k)
	{
		int delta = Integer.signum(k);
		LocatedGroup out = nbhGroup[1][1][delta+1];
		k -= delta;
		while ( k != 0 )
		{
			out = moveZ(delta);
			k -= delta;
		}
		return out;
	}
	
	/**
	 * \brief Compute distance to closest carrier.
	 */
	public void distanceFromBorders()
	{
		Double valueCarrier = Double.MAX_VALUE;
		for (AllBC aBoundary : agentGrid.domain.getAllSupportBoundaries() )
			valueCarrier = Math.min(valueCarrier, aBoundary.getDistance(cc));
		distanceFromCarrier = valueCarrier;
	}

	/**
	 * \brief Use the boundary conditions to build 3D neighbourhood reference
	 * map.
	 * 
	 * @param shovGrid	The shoving grid used to build reference map.
	 */
	protected void testNbh_3D(LocatedGroup[] shovGrid)
	{
		AllBC aBC;
		int index;
		DiscreteVector nbhDC = new DiscreteVector();
		for (int i = 0; i<3; i++) 
			for (int j = 0; j<3; j++)
				for (int k = 0; k<3; k++)
					try
					{
						// Get your supposed neighbour
						nbhDC.set(dc.i+i-1, dc.j+j-1, dc.k+k-1);
						index = agentGrid.getIndexedPosition(nbhDC);

						// If the neighbor is outside the domain, apply appropriate boundary
						// conditions until no more boundaries are crossed; the neighbor will
						// then be within the domain if appropriate (periodic boundary) or will
						// be the outside-domain neighbor that was picked originally
						if ( shovGrid[index].isOutside )
						{
							int oldindex;
							do {
								oldindex = index;
								aBC = agentGrid.domain.testCrossedBoundary(shovGrid[index].cc);
								if (aBC == null)
									break; // no boundary was crossed
								index = agentGrid.getIndexedPosition(aBC.lookAt(shovGrid[index].cc));
							} while (oldindex != index);
						}
						/*
						 * Store the reference to your neighbour.
						 */
						nbhIndex[i][j][k] = index;
						nbhGroup[i][j][k] = shovGrid[index];
					}
					catch (Exception e)
					{
						// nothing done here
					}
	}

	/**
	 * \brief Use the boundary conditions to build 3D neighbourhood reference
	 * map.
	 * 
	 * @param shovGrid	The shoving grid used to build reference map.
	 */
	protected void testNbh_2D(LocatedGroup[] shovGrid)
	{
		AllBC aBC;
		int index, oldIndex;
		DiscreteVector nbhDC = new DiscreteVector();
		ContinuousVector cc;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				try
				{
					// Get your supposed neighbour.
					nbhDC.set(dc.i + i - 1, dc.j + j - 1, dc.k);
					index = agentGrid.getIndexedPosition(nbhDC);
					// If the neighbor is outside the domain, apply appropriate boundary
					// conditions until no more boundaries are crossed; the neighbor will
					// then be within the domain if appropriate (periodic boundary) or will
					// be the outside-domain neighbor that was picked originally
					if ( shovGrid[index].isOutside )
						do {
							oldIndex = index;
							cc = shovGrid[index].cc;
							aBC = agentGrid.domain.testCrossedBoundary(cc);
							if ( aBC == null )
								break; // no boundary was crossed
							index = agentGrid.getIndexedPosition(aBC.lookAt(cc));
						} while ( oldIndex != index );
					// Store the reference to your neighbour
					nbhIndex[i][j][1] = index;
					nbhGroup[i][j][1] = shovGrid[index];
				}
				catch (Exception e)
				{
					// Nothing done here.
				}
	}

	/**
	 * \brief Calculate the number of neighbours that contains no biomass.
	 * 
	 * In 2D, we count y-side neighbours twice to be consistent with 3D.
	 * 
	 * @return the number of empty neighbours.
	 */
	public int freeNbh()
	{
		nFreeNbh = 0;
		// x-side neighbours
		if ( nbhGroup[0][1][1].status == 2 )
			nFreeNbh++;
		if ( nbhGroup[2][1][1].status == 2 )
			nFreeNbh++;
		// y-side neighbours:
		if ( nbhGroup[1][0][1].status == 2 )
			nFreeNbh++;
		if ( nbhGroup[1][2][1].status == 2 )
			nFreeNbh++;
		// z-side neighbours:
		if ( agentGrid.is3D )
		{
			if ( nbhGroup[1][1][0].status == 2 )
				nFreeNbh++;
			if ( nbhGroup[1][1][2].status == 2 )
				nFreeNbh++;
		}
		else
		{
			if ( nbhGroup[1][0][1].status == 2 )
				nFreeNbh++;
			if ( nbhGroup[1][2][1].status == 2 )
				nFreeNbh++;
		}
		return nFreeNbh;
	}
	
	/**
	 * Used for debugging purposes.
	 */
	public int countNullNeighbours()
	{
		int numNull = 0;
		for ( LocatedGroup[][] nbhZ : nbhGroup )
			for ( LocatedGroup[] nbhY : nbhZ )
				for ( LocatedGroup nbhX : nbhY )
					if ( nbhX == null )
						numNull++;
		return numNull;
	}

	/**
	 * \brief Comparator used by the detachment levelset algorithm
	 * 
	 * @author Chaodong Zhang
	 * @author Robert Clegg
	 */
	public static class TValueComparator implements java.util.Comparator<Object> 
	{
		@Override
		public int compare(Object b1, Object b2) 
		{
			Double temp = ((LocatedGroup) b1).erosionTime -
						  ((LocatedGroup) b2).erosionTime;
			temp = Math.signum(temp);
			return temp.intValue();
		}
	}
}
