/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation and their participation in relevant reactions
 * 
 * Package of utilities that create and manage agents in the simulation and their participation in relevant reactions. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;

import java.util.*;


import simulator.Simulator;
import simulator.AgentContainer;
import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.SoluteGrid;
import utils.ExtraMath;

/**
 * \brief Object to hold a group of agents in one location on the agent grid
 * 
 * Object to hold a group of agents in one location on the agent grid
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class LocatedGroup {

	/**
	 * Agent container this located group is within
	 */
	public AgentContainer           agentGrid;
	
	/**
	 * Linked list to hold members of this group
	 */
	public LinkedList<LocatedAgent> group              = new LinkedList<LocatedAgent>();

	/**
	 * Concentration of species in this group, thus area represented
	 */
	public double[]                 speciesConcentration;
	
	/**
	 * Total volume of agents in this group
	 */
	public double                   totalVolume        = 0;
	
	public double                   deltaV             = 0;
	
	/**
	 * Total concentration of agents in this group
	 */
	public double                   totalConcentration = 0;
	
	/**
	 * Total mass of agents in this group
	 */
	public double                   totalMass          = 0;
	
	public double                   erosionTime        = Double.NaN;
	
	/**
	 * Distance of the group at this location from the carrier
	 */
	public double                   distanceFromCarrier;
	
	/**
	 * Distance of the group at this location from the bulk
	 */
	public double 					distanceFromBulk;
	
	public double 					ratio;

	/**
	 * The index of the grid at which this located group represents
	 */
	public int                      gridIndex;
	
	/**
	 * Coordinates of this location as a continuous vector
	 */
	public ContinuousVector         cc;
	
	/**
	 * Coordinates of this location as a discrete
	 */
	public DiscreteVector           dc;
	
	/**
	 * Index of the neighbours of this group
	 */
	public int[][][]                nbhIndex           = new int[3][3][3];
	
	/**
	 * Neighbouring groups around this group
	 */
	public LocatedGroup[][][]       nbhGroup           = new LocatedGroup[3][3][3];

	/**
	 * Space occupation (-1->outside, 0->carrier, 1->biofilm, 2->liquid, 3->bulk)
	 */ 
	public int                      status             = 2;
	
	/**
	 * Boolean stating whether this location is in the bulk
	 */
	public boolean                  isBulk             = false;
	
	/**
	 * Boolean stating whether this location is in the carrier
	 */
	public boolean                  isCarrier          = false;
	
	/**
	 * Boolean stating whether this location is outside the grid
	 */
	public boolean                  isOutside;

	/**
	 * Number of free neighbours around this location
	 */
	public int                      nFreeNbh;
	
	/**
	 * Vector to hold an amount of distance an agent is to move
	 */
	public ContinuousVector         move               = new ContinuousVector();

	

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
		agentGrid = anAgentGrid;
		// Spatial location of the group
		gridIndex = index;
		
		ratio = 0;

		// Coordinates if padding is removed
		dc = agentGrid.getGridPosition(gridIndex);
		cc = agentGrid.getGridLocation(gridIndex);

		// Initialise biomass statistics
		speciesConcentration = new double[aSimulator.speciesDic.size()];

		// Check if the group is inside the domain
		isOutside = false;
		for (AllBC aBC : anAgentGrid.domain.getAllBoundaries()) {
			if (aBC.isOutside(cc)) {
				isOutside = true;
				status = -1;
				// setBoundary sets the status to 3!
				aBC.setBoundary(this);
				break;
			}
		}
	}

	/**
	 * \brief Builds neighbourhood reference map
	 * 
	 * Builds neighbourhood reference map
	 */
	public void init() {
		if (isOutside) return;

		if (agentGrid.is3D) testNbh_3D(agentGrid.getShovingGrid());
		else testNbh_2D(agentGrid.getShovingGrid());

		distanceFromBorders();
		if (distanceFromCarrier<agentGrid.getResolution()) isCarrier = true;
		if (distanceFromBulk<agentGrid.getResolution()) isBulk = true;
	}

	/**
	 * \brief Refresh status and concentration of the group
	 * 
	 * Refresh status and concentration of the group
	 * 
	 */
	public void refreshElement() {
		double volume, value;
		volume = ExtraMath.cube(agentGrid.getResolution());
		//sonia:chemostat
		if(Simulator.isChemostat){
			//do not refresh groups
		}else{
			// Refresh group status (carrier, biofilm, free)
			if (status>0) status = (group.size()>0 ? 1 : 2);
			if (isCarrier) status = 0;
		}
		// Refresh biomass density
		totalConcentration = 0;
		totalMass = 0;
		Arrays.fill(speciesConcentration, 0);

		for (LocatedAgent aLoc : group) {
			totalMass += aLoc.getTotalMass();
			value = aLoc.getTotalMass()/volume;

			// we treat agents as cylinders in 2D, and so the element volume is
			// the same cube as in 3D; no concentration correction is necessary

			totalConcentration += value;
			speciesConcentration[aLoc.speciesIndex] += value;
		}
	}

	/**
	 * \brief Refresh the volume statistics of this group
	 * 
	 * Refresh the volume statistics of this group
	 * 
	 * @return	Double noting the total volume of agents in this group
	 */
	public double refreshVolume() {
		totalVolume = 0;
		for (LocatedAgent aLoc : group) {
			totalVolume += aLoc.getVolume(true);
		}
		return totalVolume;
	}

	/**
	 * \brief Compute the gradient due to pressure and use it to set the advective affect
	 * 
	 * Compute the gradient due to pressure and use it to set the advective affect
	 * 
	 * @param pressure	Pressure grid
	 * @param deltaT	DeltaT
	 * @return	Norm of the movement vector under the affect of the pressure gradient
	 */
	public double computeMove(SoluteGrid pressure, double deltaT) 
	{
		if (this.isOutside) move.reset();
		else {
			move = pressure.getGradient(this.cc);
			if (!move.isValid()) move.reset();
			move.times(-deltaT);
		}
		return move.norm();
	}

	public void resetMove() {
		move.reset();
	}

	/**
	 * \brief Scale the movement vector for the grid element and apply to each agent
	 * 
	 * Scale the movement vector for the grid element and apply to each agent
	 * 
	 * @param alpha	Scaling factor to be applied to a move
	 */
	public void addMoveToAgents(double alpha) {
		for (LocatedAgent aLoc : group) {
			//if (aLoc.isAttached()) move.x = 0;
			move.times(alpha);
			aLoc.addMovement(move);
		}
	}

	/**
	 * \brief Safely remove located agent from agentList & agentGrid
	 * 
	 * Safely remove located agent from agentList & agentGrid
	 * 
	 */
	public void killAll() {
		ListIterator<LocatedAgent> iter = group.listIterator();
		LocatedAgent aLoc;
		while (iter.hasNext()) {
			aLoc = iter.next();
			//sonia 27.04.2010
			//added add aLoc to agentToKill list and reason of death;
			aLoc.death="detachment";
			agentGrid._agentToKill.add(aLoc);
			agentGrid.agentList.remove(aLoc);
			iter.remove();
		}
	}

	/**
	 * \brief Remove an agent from this LocatedGroup
	 * 
	 * Remove an agent from this LocatedGroup
	 * @param anAgent	LocatedAgent to remove from this group
	 */
	public void remove(LocatedAgent anAgent) {
		group.remove(anAgent);
		//sonia:chemostat

		if(Simulator.isChemostat){

		}else{
			if (group.isEmpty()) status = 2;
		}
	}

	/**
	 * \brief Add an agent to this LocatedGroup
	 * 
	 * Add an agent to this LocatedGroup
	 * @param anAgent	LocatedAgent to add to this group
	 */
	public void add(LocatedAgent anAgent) 
	{
		group.add(anAgent);
		status = 1;
		anAgent.setGridIndex(gridIndex);
	}

	/**
	 * \brief Shuffle coordinates of an agent inside the group
	 * 
	 * Shuffle coordinates of an agent inside the group.
	 * 
	 * @param aLoc	LocatedAgent which is to be shuffled
	 */
	public void host(LocatedAgent aLoc) {
		double res = agentGrid.getResolution();
		ContinuousVector cc = aLoc.getLocation();
		cc.y = ExtraMath.getUniRand(this.cc.y-res/2, this.cc.y+res/2);
		cc.z = ExtraMath.getUniRand(this.cc.z-res/2, this.cc.z+res/2);
		cc.x = this.cc.x-res/2;
		add(aLoc);
	}

	/**
	 * \brief Move the X parameter by a specified amount
	 * 
	 * Move the X parameter by a specified amount
	 * 
	 * @param i	The current I grid element
	 * @return	Located group reached by that move
	 */
	public LocatedGroup moveX(int i) {
		int delta = (int) Math.signum(i);
		LocatedGroup out = nbhGroup[delta+1][1][1];
		i -= delta;
		while (i!=0) {
			out = moveX(delta);
			i -= delta;
		}
		return out;
	}

	/**
	 * \brief Move the Y parameter by a specified amount
	 * 
	 * Move the Y parameter by a specified amount
	 * 
	 * @param i	The current grid element
	 * @return	Located group reached by that move
	 */
	public LocatedGroup moveY(int i) {
		int delta = (int) Math.signum(i);
		LocatedGroup out = nbhGroup[1][delta+1][1];
		i -= delta;
		while (i!=0) {
			out = moveY(delta);
			i -= delta;
		}
		return out;
	}

	/**
	 * \brief Move the Z parameter by a specified amount
	 * 
	 * Move the Z parameter by a specified amount
	 * 
	 * @param i	The current grid element
	 * @return	Located group reached by that move
	 */
	public LocatedGroup moveZ(int i) {
		int delta = (int) Math.signum(i);
		LocatedGroup out = nbhGroup[1][1][delta+1];
		i -= delta;
		while (i!=0) {
			out = moveZ(delta);
			i -= delta;
		}
		return out;
	}

	/**
	 * \brief Compute the difference vector between two continuous vectors
	 * 
	 * Compute the difference vector between two continuous vectors
	 * 
	 * @param me	Location of one agent expressed as continuous vector
	 * @param him	Location of second agent expressed as continuous vector
	 * @param move	ContinuousVector to hold the difference between the two vectors
	 * @return	Double distance between the two points
	 * @deprecated
	 */
	public double computeDifferenceVector(ContinuousVector me, ContinuousVector him,
			ContinuousVector move) {
		double gridLength;

		move.x = me.x-him.x;
		// check periodicity in X
		gridLength = agentGrid.domain.length_X;
		if (Math.abs(move.x)>.5*gridLength) {
			move.x -= Math.signum(move.x)*gridLength;
		}

		move.y = me.y-him.y;
		// check periodicity in Y
		gridLength = agentGrid.domain.length_Y;
		if (Math.abs(move.y)>.5*gridLength) {
			move.y -= Math.signum(move.y)*gridLength;
		}

		if (agentGrid.is3D) {
			move.z = me.z-him.z;
			// check periodicity in Z
			gridLength = agentGrid.domain.length_Z;
			if (Math.abs(move.z)>.5*gridLength) {
				move.z -= Math.signum(move.z)*gridLength;
			}

		} else {
			move.z = 0;
		}
		double d = Math.sqrt(move.x*move.x+move.y*move.y+move.z*move.z);

		return d;
	}

	/**
	 * \brief Compute distance to closest carrier and closest bulk
	 * 
	 * Compute distance to closest carrier and closest bulk
	 */
	public void distanceFromBorders() {
		LinkedList<AllBC> allBoundary = agentGrid.domain.getAllBoundaries();

		double valueCarrier = Double.MAX_VALUE;
		double valueBulk = Double.MAX_VALUE;

		for (AllBC aBoundary : allBoundary) {
			if (aBoundary.isSupport()) valueCarrier = Math.min(valueCarrier, aBoundary
					.getDistance(cc));
			if (aBoundary.hasBulk()) valueBulk = Math.min(valueBulk, aBoundary.getDistance(cc));
		}
		distanceFromCarrier = valueCarrier;
		distanceFromBulk = valueBulk;
	}

	/**
	 * \brief Use the boundary conditions to build 3D neighbourhood reference map
	 * 
	 * Use the boundary conditions to build 3D neighbourhood reference map
	 * 
	 * @param shovGrid	The shoving grid used to build reference map
	 */
	protected void testNbh_3D(LocatedGroup[] shovGrid) {
		AllBC aBC;
		int index;
		DiscreteVector nbhDC = new DiscreteVector();

		// Build neighbourhood reference map
		for (int i = 0; i<3; i++) {
			for (int j = 0; j<3; j++) {
				for (int k = 0; k<3; k++) {
					try {
						// Get your supposed neighbour
						nbhDC.set(dc.i+i-1, dc.j+j-1, dc.k+k-1);
						index = agentGrid.getIndexedPosition(nbhDC);

						// If the neighbor is outside the domain, apply appropriate boundary
						// conditions until no more boundaries are crossed; the neighbor will
						// then be within the domain if appropriate (periodic boundary) or will
						// be the outside-domain neighbor that was picked originally
						if (shovGrid[index].isOutside) {
							int oldindex;
							do {
								oldindex = index;
								aBC = agentGrid.domain.testCrossedBoundary(shovGrid[index].cc);
								if (aBC == null) break; // no boundary was crossed
								index = agentGrid.getIndexedPosition(aBC.lookAt(shovGrid[index].cc));
							} while (oldindex != index);
						}

						// Store the reference to your neighbour
						nbhIndex[i][j][k] = index;
						nbhGroup[i][j][k] = shovGrid[index];
					} catch (Exception e) {
						// nothing done here

					}
				}
			}
		}
	}

	/**
	 * \brief Use the boundary conditions to build 3D neighbourhood reference map
	 * 
	 * Use the boundary conditions to build 3D neighbourhood reference map
	 * 
	 * @param shovGrid	The shoving grid used to build reference map
	 */
	protected void testNbh_2D(LocatedGroup[] shovGrid) 
	{
		AllBC aBC;
		int index;
		DiscreteVector nbhDC = new DiscreteVector();

		// Build neighbourhood reference map
		for (int i = 0; i<3; i++) {
			for (int j = 0; j<3; j++) {
				try {
					// Get your supposed neighbour
					nbhDC.set(dc.i+i-1, dc.j+j-1, dc.k);
					index = agentGrid.getIndexedPosition(nbhDC);

					// If the neighbor is outside the domain, apply appropriate boundary
					// conditions until no more boundaries are crossed; the neighbor will
					// then be within the domain if appropriate (periodic boundary) or will
					// be the outside-domain neighbor that was picked originally
					if (shovGrid[index].isOutside) {
						int oldindex;
						do {
							oldindex = index;
							aBC = agentGrid.domain.testCrossedBoundary(shovGrid[index].cc);
							if (aBC == null) break; // no boundary was crossed
							index = agentGrid.getIndexedPosition(aBC.lookAt(shovGrid[index].cc));
						} while (oldindex != index);
					}

					// Store the reference to your neighbour
					nbhIndex[i][j][1] = index;
					nbhGroup[i][j][1] = shovGrid[index];
				} catch (Exception e) {
					// nothing done here

				}
			}
		}
	}

	/**
	 * \brief Calculate the number of neighbours that contains no biomass
	 * 
	 * Calculate the number of neighbours that contains no biomass
	 * 
	 * @return the number of empty neighbours
	 */
	public int freeNbh() {
		nFreeNbh = 0;
		// Rob Feb 2011: add y-side neighbours twice in 2D
		// (to be consistent with 3D)

		// x neighbours
		if (nbhGroup[0][1][1].status==2) nFreeNbh++;
		if (nbhGroup[2][1][1].status==2) nFreeNbh++;

		// y-side neighbours:
		if (nbhGroup[1][0][1].status==2) nFreeNbh++;
		if (nbhGroup[1][2][1].status==2) nFreeNbh++;

		// z-side neighbours:
		if (agentGrid.is3D) {
			if (nbhGroup[1][1][0].status==2) nFreeNbh++;
			if (nbhGroup[1][1][2].status==2) nFreeNbh++;
		} else {
			if (nbhGroup[1][0][1].status==2) nFreeNbh++;
			if (nbhGroup[1][2][1].status==2) nFreeNbh++;
		}
		return nFreeNbh;
	}

	/**
	 * \brief Comparator used by the detachment levelset algorithm
	 * 
	 * Comparator used by the detachment levelset algorithm
	 * 
	 * @author Chaodong Zhang
	 */
	public static class TValueComparator implements java.util.Comparator<Object> 
	{

		public int compare(Object b1, Object b2) 
		{
			double f1=((LocatedGroup) b1).erosionTime;
			double f2=((LocatedGroup) b2).erosionTime;
			
			if(f1==f2)
				return 0;
			else if(f1 > f2)
				return 1;
			else
				return -1;
		}
	}

	/**
	 * \brief Comparator used by the shrinking levelset algorithm
	 * 
	 * Comparator used by the shrinking levelset algorithm
	 * 
	 * @author lal
	 */
	public static class DistanceValueComparator implements java.util.Comparator<LocatedGroup> {

		public int compare(LocatedGroup b1, LocatedGroup b2) {
			double out = b1.distanceFromCarrier-b2.distanceFromCarrier;

			if (out==0) {
				return (b1.deltaV>b2.deltaV ? 1 : -1);
			} else if (out>0) {
				return 1;
			} else {
				return -1;
			}

		}
	}
}
