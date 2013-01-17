/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  ___________________________________________________________________________
 * IsBoundaryCondition : Interface for the boundary conditions on the system's 
 * margins
 * To be used in solvers and agent move
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.geometry.boundaryConditions;

import org.jdom.Element;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.geometry.*;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;

public interface IsBC {

	/* Functions related to the geometry of the boundary __________________ */

	/**
     * Used during the initialisation, load the class describing the shape of
     * the boundary defined in the parent class
     * @param Element
     */
	public void readGeometry(Element geometryRoot, Domain aDomain);

	// Geometric methods, defined in the parent class
	public double getDistance(ContinuousVector cc);

	public ContinuousVector getIntersection(ContinuousVector position, ContinuousVector vector);

	public ContinuousVector getOrthoProj(ContinuousVector cc);

	public boolean isOutside(ContinuousVector cc);
	public boolean isOutside(DiscreteVector dC, SpatialGrid aSpatialGrid);
	
	public boolean isActive();
	public boolean isCyclic();

	public boolean hasBulk();

	public double updateBulk(SoluteGrid[] reacGrid,double time,double timeStepMin,boolean lastIter);
	public double getBulkValue(int soluteIndex);
	public Bulk getBulk();
	public boolean isSupport();

	/**
     * Initialisation of the boundary
     */
	public void init(Simulator aSim, Domain aDomain, Element aBoundCondMarkUp);

	/* __________________ INTERACTION WITH THE SOLVER _____________________ */

	/**
     * Refreshes the boundary conditions during multigrid computatioin
     * @see ComputationDomain.refreshBoundaries()
     */
	public void refreshBoundary(SoluteGrid aSoluteGrid);
	public void refreshDiffBoundary(SoluteGrid aSoluteGrid,SoluteGrid aSolutegrid);

	/* ______________INTERACTION WITH THE PARTICLES _____________________ */

	public ContinuousVector lookAt(ContinuousVector cc);
	public void setBoundary(LocatedGroup aGroup);
	
	/**
     * Modify the movement vector : the new position is the orthognal projection
     * on the boundary surface
     * @see LocatedAgent.move();
     */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc);

}