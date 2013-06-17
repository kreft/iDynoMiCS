/**
 * \package simulator.geometry.boundaryConditions
 * \brief Package of boundary conditions that can be used to capture agent behaviour at the boundary of the computation domain
 * 
 * Package of boundary conditions that can be used to capture agent behaviour at the boundary of the computation domain. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry.boundaryConditions;

import utils.XMLParser;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;

/**
 * \brief BoundaryConstant : the concentration on the boundary is fixed by a constant bulk, the agents crossing this line die
 * 
 * BoundaryConstant : the concentration on the boundary is fixed by a constant bulk, the agents crossing this line die. The boundary 
 * represents, for example, the connection to a larger system where the concentration can be considered constant. The solute 
 * concentrations at this boundary are fixed, and agents crossing this boundary are considered to have entered the planktonic bulk 
 * domain. 
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class BoundaryConstant extends AllBC{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The reference to the bulk which maintains the concentration
	 */
	private Bulk              _connectedBulk;

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file, and loads the description of the connected bulk
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file, and loads the description of the connected bulk
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBoundCondMarkUp	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) {
		// Load the geometry of the boundary
		readGeometry(aBoundCondMarkUp, aDomain);
		aDomain.addBoundary(this);
		
		// Load description of the connected bulk
		_connectedBulk = aSim.world.getBulk(aBoundCondMarkUp.getParam("bulk"));		
	}

	/**
	 * \brief Solver for the constant concentration boundary condition. Initialises the course along the shape of the boundary, setting the values of solute near the boundary as required 
	 * 
	 * Solver for the variable concentration boundary condition. Initialises the course along the shape of the boundary, , setting the values of solute near the boundary as required
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) {
		// Some internal variables
		double bulkValue = _connectedBulk.getValue(aSoluteGrid.soluteIndex);

		// Initialise the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
			aSoluteGrid.setValueAt(bulkValue, dcOut);
		}
	}

	/**
	 * \brief Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions need the input corrected, some don't and just return the input
	 * 
	 * Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions (e.g. BoundaryCyclic_ 
	 * need the input corrected due to the condition, some don't and just return the input. Maybe we'll change this at some point as to 
	 * just return the input looks a bit daft - but we'll leave it here for the moment
	 * 
	 * @param cc	ContinuousVector that gives the current location of an agent to check on the grid
	 */
	public ContinuousVector lookAt(ContinuousVector cc) {
		return cc;
	}

	/**
     * \brief Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * @param aGroup	LocatedGroup object which has been detected to be outside the boundary
     */
	public void setBoundary(LocatedGroup aGroup) {
		aGroup.status = 3;
		// status 3 -> bulk
	}

	/**
	 * \brief Kills any agents that are crossing this boundary as they are leaving the simulated system
	 * 
	 * Kills any agents that are crossing this boundary as they are leaving the simulated system
	 * 
	 * @param anAgent	The LocatedAgent that has crossed the boundary
	 * @param target	Vector of where this agent was going to be placed
	 */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) {
		
		//sonia 27.04.2010
		//recording reason of death (agent will be moved to agentToKill list when die() calls registerDeath()
		anAgent.death = "overBoard";
		
		anAgent.die(false);
		// to label this agent as "shoving solved", set to zero its movement.
		anAgent.getMovement().reset();
		target.set(anAgent.getLocation());
	}

}
