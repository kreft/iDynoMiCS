/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 * ___________________________________________________________________________
 * BoundaryBulk : the concentration on the boundary is fixed by a dynamical 
 * bulk, the agents crossing this line die
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */

package simulator.geometry.boundaryConditions;

import utils.XMLParser;

import simulator.*;
import simulator.geometry.*;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;

public class BoundaryBulk extends AllBC{

	/* _____________________________ FIELDS _______________________________ */
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	// At which bulk the liquid phase is connected
	protected Bulk            _connectedBulk;
	static double             bulkValue;

	/* ___________________________ CONSTRUCTOR _____________________________ */
	public BoundaryBulk() {
		hasBulk = true;
	}

	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) {
		// Load the geometry of the boundary
		readGeometry(aBoundCondMarkUp, aDomain);		
		
		aDomain.addBoundary(this);
		
		// Load description of the connected bulk
		String bulkName = aBoundCondMarkUp.getParam("bulk");
		_connectedBulk = aSim.world.getBulk(bulkName);
	}

	/* ____________________________ SOLVER _________________________________ */

	public void refreshBoundary(SoluteGrid aSoluteGrid) {
		// Store the concentration in the bulk
		bulkValue = _connectedBulk.getValue(aSoluteGrid.soluteIndex);
		// bulkValue = 1;
		// Initialise the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
			aSoluteGrid.setValueAt(bulkValue, dcOut);
		}
	}

	// bvm note 15.12.08: modified inputs to allow reaction or flux-based bulk treatment
	public void updateBulk(SoluteGrid[] allSG, SoluteGrid[] allRG, double timeStep) {
		_connectedBulk.updateBulk(allSG, allRG, timeStep);
	}

	public Bulk getBulk() {
		return _connectedBulk;
	}

	public double getBulkValue(int soluteIndex) {
		return _connectedBulk.getValue(soluteIndex);
	}

	/* _______________________ LOCATED Agents ______________________________ */

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
		//LogFile.writeLog("agent killed by Bulk Boundary");
		
		// to label this agent as "shoving solved", set to zero its movement.
		anAgent.getMovement().reset();
		target.set(anAgent.getLocation());
		
	}



	public String toString() {
		return new String("Bulk:"+this._mySide);
	}
}
