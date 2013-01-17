/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *______________________________________________________
 * AllBoundaryCondition : group all methods expected by the interface but 
 * common to most of the boundary classes
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.geometry.boundaryConditions;

import utils.XMLParser;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import simulator.geometry.shape.IsShape;

public abstract class AllBC{

	/* _______________________________ FIELDS _________________________________ */
	// The name of the boundary describing its side (xOy,...)
	protected String                _mySide;
	// The shape of the boundary
	protected IsShape               _myShape;
	// Quick test to know if the boundary is cyclic
	protected boolean               isCyclic        = false;
	protected boolean               _isSupport      = false;
	protected boolean               hasBulk         = false;
	protected boolean               activeForSolute = true;
	protected boolean               hasAgar         = false;
	
	/* ____________________ INTERNAL TEMPRARY VARIABLES ________________________ */
	// Discrete coordinates of a voxel inside the computation domain but along
	// the boundary
	protected static DiscreteVector dcIn            = new DiscreteVector();
	// Discrete coordinates of the voxel in front of the other one
	protected static DiscreteVector dcOut           = new DiscreteVector();

	/* ________________________ CONSTRUCTION METHODS ________________________ */

	/**
	 * Generic constructor called to dynamically instanciate a child class object
	 */
	public static AllBC staticBuilder(XMLParser root, Simulator aSim,Domain aDomain) {
		// Create the object
		AllBC out = (AllBC) root.instanceCreator("simulator.geometry.boundaryConditions");
		
		// Initialise & declare the boundary
		out.init(aSim, aDomain, root);

		return out;
	}
	
	public abstract void init(Simulator aSim,Domain aDomain,XMLParser root);

	
	/**
     * Used during the initialisation, load the class describing the shape of
     * the boundary
     * @param Element
     */
	public void readGeometry(XMLParser geometryRoot, Domain aDomain) {
		// Set the name of the boundary
		_mySide = geometryRoot.getAttributeStr("name");

		// Set the class to use to define the shape
		String className = "simulator.geometry.shape.";
		className += geometryRoot.getChild("shape").getAttributeStr("class");

		// Build the instance used to describe the shape
		try {
			_myShape = (IsShape) Class.forName(className).newInstance();
			_myShape.readShape(new XMLParser(geometryRoot.getChildElement("shape")), aDomain);

		} catch (Exception e) {
		}
	}
	
	
	public boolean isOutside(ContinuousVector cc) {
		return _myShape.isOutside(cc);
	}
	




	/* __________________ INTERACTION WITH THE SOLVER _____________________ */

	/**
     * Refreshes the boundary conditions during Diff&Reac solving
     * @see ComputationDomain.refreshBoundaries()
     */
	public abstract void refreshBoundary(SoluteGrid aSoluteGrid);
	
	/**
     * Method used if a boundary modifies the local diffusivity constant
     * Most of boundaries do not modify it
     * @param aSoluteGrid
     */
	public void refreshDiffBoundary(SoluteGrid aSoluteGrid,SoluteGrid aSolutegrid){};

	/* ______________INTERACTION WITH THE PARTICLES _____________________ */

	public abstract ContinuousVector lookAt(ContinuousVector cc);
	public abstract void setBoundary(LocatedGroup aGroup);
	
	/**
     * Modify the movement vector : the new position is the orthognal projection
     * on the boundary surface
     * @see LocatedAgent.move();
     */
	public abstract void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc);
	
	
	
	


	/* ___________________ INTERACTION WITH THE DOMAIN _________________________ */

	public boolean isCyclic() {
		return isCyclic;
	}

	public boolean isSupport() {
		return _isSupport;
	}

	public boolean isActive(){
		return activeForSolute;
	}
	
	public boolean hasBulk() {
		return hasBulk;
	}

	public void updateBulk(SoluteGrid[] soluteGrid, SoluteGrid[] reacGrid, double timeStep) {};

	public boolean hasAgar() {
		return hasAgar;
	}

	public void updateAgar(SoluteGrid[] soluteGrid,SoluteGrid[] reactionGrid, double timeStep) {};

	public boolean isOutside(DiscreteVector dc, SpatialGrid aSpatialGrid) {
		return _myShape.isOutside(new ContinuousVector(dc, aSpatialGrid.getResolution()));
	}

	/* ____________________ TOOLBOX ______________________________ */

	/**
     * @return the intersection between the boundary and a line (defined by a
     * position and a vector)
     */
	public ContinuousVector getIntersection(ContinuousVector position, ContinuousVector vector) {
		return _myShape.intersection(position, vector);
	}

	/**
     * @return the orthogonal projection of a location on the boundary
     * @param cc
     * @return
     */
	public ContinuousVector getOrthoProj(ContinuousVector cc) {
		return _myShape.getOrthoProj(cc);
	}

	public Bulk getBulk() {
		return null;
	}

	public IsShape getShape() {
		return _myShape;
	}
	
	public String getSide() {
		return _mySide;
	}
	
	/**
     * @param cc
     * @return the distance between a location and the boundary
     */
	public double getDistance(ContinuousVector cc) {
		return _myShape.getDistance(cc);
	}

	public double getBulkValue(int soluteIndex) {
		return 0;
	}
	

}
