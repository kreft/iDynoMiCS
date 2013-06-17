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

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import utils.XMLParser;

/**
 * \brief Defines an impermeable boundary to solutes and agents
 * 
 * Defines an impermeable boundary to solutes and agents. As a consequence the normal components of solute concentration gradients 
 * will be zero at this boundary. Agents attempting to cross the boundary are prevented from doing so
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BoundaryZeroFlux  extends AllBC
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Stores a vector normal to the boundary and starting from the orthogonal projection
	 */
	private static ContinuousVector vectorIn;

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBoundCondMarkUp	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) 
	{
		readGeometry(aBoundCondMarkUp, aDomain);
		aDomain.addBoundary(this);
		_isSupport = true;
	}

	/* ________________________ SOLVER _______________________________ */
	
	/**
	 * \brief Solver for the zero flux boundary condition. Initialises the course along the shape of the boundary and moves any points outside the domain inside
	 * 
	 * Solver for the zero flux boundary condition. Initialises the course along the shape of the boundary and moves any points outside the domain inside
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) 
	{
		
		// Initialize the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		// Send a point belonging to the boundary and the closest point outside
		// the domain
		while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) 
		{
			aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
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
	public void setBoundary(LocatedGroup aGroup) 
	{
		// status 0 -> carrier
		aGroup.status = 0;
	}
	
	/**
	 * \brief Applies the boundary condition by modifying the movement vector. New position is orthogonal projection of the outside point on the boundary surface
	 * 
	 * Applies the boundary condition by modifying the movement vector. New position is orthogonal projection of the outside point on the boundary surface
	 * 
	 * @param anAgent	The Located Agent which is attempting to cross the boundary
	 * @param target	The target position that the agent is moving to
	 */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) {
		// Define coordinates of the corrected position
		_myShape.orthoProj(target,target);
		
		// Build a vector normal to the boundary and starting from the
        // orthogonal projection		
		vectorIn = new ContinuousVector(_myShape.getNormalInside(target));
		
		// The whole cell has to be inside, so make a translation equal to the
        // total radius		
		vectorIn.times(anAgent.getRadius(true));
		
		// Compute the new position
		target.add(vectorIn);
		
		// Compute and update the movement vector leading to this new position
		anAgent.getMovement().sendDiff(anAgent.getLocation(), target);
	}
	
	/**
	 * \brief Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * @return String noting the side of the domain that this condition applies to (i.e. x0z, xNz, etc)
	 */
	public String toString(){
		return new String("ZeroFlux:"+this._mySide);
	}
	

}
