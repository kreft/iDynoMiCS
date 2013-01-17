/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 *  ___________________________________________________________________________
 * BoundaryZeroFlux : defines an impermeable boundary
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

package simulator.geometry.boundaryConditions;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import utils.XMLParser;

public class BoundaryZeroFlux  extends AllBC{
	
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	/* ____________________ INTERNAL TEMPRARY VARIABLES ________________________ */
	private static ContinuousVector vectorIn;

	
	/* ________________________ CONSTRUCTOR _______________________________ */
	public BoundaryZeroFlux() {
	}

	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) {
		readGeometry(aBoundCondMarkUp, aDomain);
		aDomain.addBoundary(this);
		_isSupport = true;
	}

	/* ________________________ SOLVER _______________________________ */
	
	public void refreshBoundary(SoluteGrid aSoluteGrid) {
		// Initialize the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		// Send a point belonging to the boundary and the closest point outside
		// the domain
		while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
			aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
		}
	}
	

	/* _______________________ FOR LOCATED AGENTS __________________________ */
	/**
	 * 
	 */
	public ContinuousVector lookAt(ContinuousVector cc) {
		return cc;
	}

	/**
     * Label a LocatedGroup which has been identified being outside this boundary
     */
	public void setBoundary(LocatedGroup aGroup) {
		// status 0 -> carrier
		aGroup.status = 0;
	}
	
	/**
     * Modify the movement vector : the new position is the orthognal projection
     * of the outside point on the boundary surface
     * @see LocatedAgent.move();
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
	
	public String toString(){
		return new String("ZeroFlux:"+this._mySide);
	}
	

}
