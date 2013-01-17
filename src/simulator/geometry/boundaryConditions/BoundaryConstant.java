/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *_________________________________________________________________________
 * BoundaryConstant : the concentration on the boundary is fixed by a constant 
 * bulk, the agents crossing this line die
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */

package simulator.geometry.boundaryConditions;

import utils.XMLParser;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;

public class BoundaryConstant extends AllBC{

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	// The reference to the bulk which maintains the concentration
	private Bulk              _connectedBulk;

	/* ________________________ CONSTRUCTOR _______________________________ */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) {
		// Load the geometry of the boundary
		readGeometry(aBoundCondMarkUp, aDomain);
		aDomain.addBoundary(this);
		
		// Load description of the connected bulk
		_connectedBulk = aSim.world.getBulk(aBoundCondMarkUp.getParam("bulk"));		
	}

	/* _________________________ SOLVER __________________________________ */

	/**
	 * Apply the boundary condition and modifies the padded grid cells
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

	/* _____________________________ AGENTS ________________________________ */
	/**
	 * 
	 */
	public ContinuousVector lookAt(ContinuousVector cc) {
		return cc;
	}

	/**
	 * Label a LocatedGroup which has been identified being outside this
	 * boundary
	 */
	public void setBoundary(LocatedGroup aGroup) {
		aGroup.status = 3;
		// status 3 -> bulk
	}

	/**
	 * An agent is crossing the boundary ; he is leaving the simulated system,
	 * kill him
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
