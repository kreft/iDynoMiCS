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
import simulator.SpatialGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import simulator.geometry.shape.IsShape;

/**
 * \brief Group all methods expected by the interface but common to most of the boundary classes
 * 
 * Group all methods expected by the interface but common to most of the boundary classes
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public abstract class AllBC{

	/* _______________________________ FIELDS _________________________________ */
	/**
	 * The name of the boundary describing its side (xOy,...)
	 */
	protected String                _mySide;
	
	/**
	 * The shape of the boundary
	 */
	protected IsShape               _myShape;
	
	/**
	 * Quick test to know if the boundary is cyclic
	 */
	protected boolean               isCyclic        = false;
	
	/**
	 * Boolean noting whether this boundary is the supporting structure (substratum)
	 */
	protected boolean               _isSupport      = false;
	
	/**
	 * Boolean noting whether this boundary has an attached bulk
	 */
	protected boolean               hasBulk         = false;
	
	/**
	 * Boolean noting whether this boundary can contain active solute
	 */
	protected boolean               activeForSolute = true;
	
	/**
	 * Boolean noting whether this boundary is modelling agar
	 */
	protected boolean               hasAgar         = false;
	
	/* ____________________ INTERNAL TEMPRARY VARIABLES ________________________ */
	/**
	 * Discrete coordinates of a voxel inside the computation domain but along the boundary
	 */ 
	protected static DiscreteVector dcIn            = new DiscreteVector();
	
	/**
	 * Discrete coordinates of the voxel in front of the one outside the boundary
	 */
	protected static DiscreteVector dcOut           = new DiscreteVector();

	/* ________________________ CONSTRUCTION METHODS ________________________ */

	/**
	 * \brief Generic constructor called to dynamically instantiate a child class object
	 * 
	 * Generic constructor called to dynamically instantiate a child class object
	 * 
	 * @param root	Set of XML tags relating to one boundary condition
	 * @param aSim	The current simulation object used to simulate the conditions specified in this protocol file
	 * @param aDomain	The computation domain to which this boundary condition is assigned
	 */
	public static AllBC staticBuilder(XMLParser root, Simulator aSim,Domain aDomain) 
	{
		// Create the object
		AllBC out = (AllBC) root.instanceCreator("simulator.geometry.boundaryConditions");
		
		// Initialise & declare the boundary
		out.init(aSim, aDomain, root);

		return out;
	}
	
	/**
	 * \brief Initialises the boundary condition. Should be overriden by each boundary condition class file
	 * 
	 * Initialises the boundary condition. This method should be overriden by each boundary condition class file
	 * 
	 * @param aSim	The current simulation object used to simulate the conditions specified in this protocol file
	 * @param aDomain	The computation domain to which this boundary condition is assigned
	 * @param root	Set of XML tags relating to one boundary condition
	 */
	public abstract void init(Simulator aSim,Domain aDomain,XMLParser root);

	
	/**
     * \brief Used during the initialisation, load the class describing the shape of the boundary defined in the parent class
     * 
     * Used during the initialisation, load the class describing the shape of the boundary defined in the parent class
     * 
     * @param geometryRoot	Usually an XML set of elements that describe the boundary to be created
     * @param aDomain	The computational domain which this boundary is associated with
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
	
	/**
	 * \brief Determines if a point is outside the boundary
	 * 
	 * Determines if a point is outside the boundary
	 * 
	 * @param cc	ContinuousVector to check
	 * @return	Boolean value noting whether this coordinate is outside the boundary (true) or not (false)
	 */
	public boolean isOutside(ContinuousVector cc) {
		return _myShape.isOutside(cc);
	}
	
	/**
	 * \brief Return the name of the side of the domain which this boundary is on
	 * 
	 * Return the name of the side of the domain which this boundary is on
	 * 
	 * @return	String containing the name of the side of the domain (e.g. x0z, xNz, etc)
	 */
	public String getSideName()
	{
		return this._mySide;
	}


	/**
	 * \brief Solver for then boundary condition. Initialises the course along the shape of the boundary during multigrid computation 
	 * 
	 * Solver for the variable concentration boundary condition. Initialises the course along the shape of the boundary during multigrid computation
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 * see ComputationDomain.refreshBoundaries()
	 */
	public abstract void refreshBoundary(SoluteGrid aSoluteGrid);
	
	/**
     * \brief Method used if a boundary modifies the local diffusivity constant. Most of boundaries do not modify it
     * 
     * Method used if a boundary modifies the local diffusivity constant. Most of boundaries do not modify it
     * 
     * @param relDiff	Relative difference grid
     * @param aSolutegrid	Grid of solute information which is to be refreshed by the solver
     * 
     */
	public void refreshDiffBoundary(SoluteGrid relDiff,SoluteGrid aSolutegrid){};

	/* ______________INTERACTION WITH THE PARTICLES _____________________ */

	/**
	 * \brief Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions need the input corrected, some don't and just return the input
	 * 
	 * Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions (e.g. BoundaryCyclic_ 
	 * need the input corrected due to the condition, some don't and just return the input. Maybe we'll change this at some point as to 
	 * just return the input looks a bit daft - but we'll leave it here for the moment
	 * 
	 * @param cc	ContinuousVector that gives the current location of an agent to check on the grid
	 */
	public abstract ContinuousVector lookAt(ContinuousVector cc);
	
	/**
     * \brief Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * @param aGroup	LocatedGroup object which has been detected to be outside the boundary
     */
	public abstract void setBoundary(LocatedGroup aGroup);
	
	/**
     * Modify the movement vector : the new position is the orthognal projection
     * on the boundary surface
     * @see LocatedAgent.move();
     */
	public abstract void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc);


	/* ___________________ INTERACTION WITH THE DOMAIN _________________________ */

	/**
	 * \brief Determine if this boundary is cyclic
	 * 
	 * Determine if this boundary is cyclic
	 * 
	 * @return	Boolean noting whether this boundary is cyclic (true) or not (false)
	 */
	public boolean isCyclic() {
		return isCyclic;
	}

	/**
	 * \brief Determine whether this boundary is the supporting structure (substratum)
	 * 
	 * Determine whether this boundary is the supporting structure (substratum)
	 * 
	 * @return Boolean noting whether this boundary is the supporting structure (true) or not (false)
	 */
	public boolean isSupport() {
		return _isSupport;
	}

	/**
	 * \brief Determine if this boundary is active for solute
	 * 
	 * Determine if this boundary is active for solute
	 * 
	 * @return	Boolean noting whether this boundary is active for solute (true) or not (false)
	 */
	public boolean isActive(){
		return activeForSolute;
	}
	
	/**
	 * \brief Determine if this boundary is attached to a bulk
	 * 
	 * Determine if this boundary  is attached to a bulk
	 * 
	 * @return	Boolean noting whether this boundary  is attached to a bulk (true) or not (false)
	 */
	public boolean hasBulk() {
		return hasBulk;
	}

	/**
	 * \brief Updates the levels in the bulk. Allows reaction or flux-based bulk treatment
	 * 
	 *  Updates the levels in the bulk. Allows reaction or flux-based bulk treatment (BVM 151208)
	 * 
	 *  @param soluteGrid	Array of all solute grids
	 *  @param reacGrid	Array of all reaction grids	
	 *  @param timeStep	The internal timestep currently being applied in this simulation
	 *  
	 */
	public void updateBulk(SoluteGrid[] soluteGrid, SoluteGrid[] reacGrid, double timeStep) {};

	/**
	 * \brief Determine if the boundary condition is modelling agar
	 * 
	 * Determine if the boundary condition is modelling agar
	 * 
	 * @return	Boolean noting whether this condition is modelling agar (true) or not (false)
	 */
	public boolean hasAgar() {
		return hasAgar;
	}

	/**
	 * \brief If modelling an agar plate boundary, this method updates the boundary
	 * 
	 * @param soluteGrid	Grid of all solutes
	 * @param reactionGrid	Grid of all reactions
	 * @param timeStep	Current internal timestep of the simulation
	 */
	public void updateAgar(SoluteGrid[] soluteGrid,SoluteGrid[] reactionGrid, double timeStep) {};

	/**
	 * \brief Determines if a discrete vector location is outside the boundary
	 * 
	 * Determines if a discrete vector location is outside the boundary
	 * 
	 * @param dc	DiscreteVector to check
	 * @param aSpatialGrid	The grid to check whether a point is outside
	 * @return	Boolean value noting whether this coordinate is outside the boundary (true) or not (false)
	 */
	public boolean isOutside(DiscreteVector dc, SpatialGrid aSpatialGrid) {
		return _myShape.isOutside(new ContinuousVector(dc, aSpatialGrid.getResolution()));
	}

	/* ____________________ TOOLBOX ______________________________ */

	/**
     * \brief Calculate the intersection between the boundary and a line (defined by a position and a vector)
     * 
     * Calculate the intersection between the boundary and a line (defined by a position and a vector)
     * 
     * @param position	A continuous vector stating the point to be used in the calculation
     * @param vector	A continuous vector stating the line to be used in the calculation	
     * @return ContinuousVector stating the point of intersection between the boundary and a line
     */
	public ContinuousVector getIntersection(ContinuousVector position, ContinuousVector vector) {
		return _myShape.intersection(position, vector);
	}

	/**
     * \brief Calculate the orthogonal projection of a location on the boundary
     * 
     * Calculate the orthogonal projection of a location on the boundary
     * 
     * @param cc	A continuous vector stating the point to be used in the calculation
     * @return ContinuousVector stating the point on the boundary after the orthogonal projection 
     */
	public ContinuousVector getOrthoProj(ContinuousVector cc) {
		return _myShape.getOrthoProj(cc);
	}

	/**
	 * \brief Return the bulk that is connected to this boundary
	 * 
	 * Return the bulk that is connected to this boundary
	 * 
	 * @return Bulk object that is connected to this boundary
	 */
	public Bulk getBulk() {
		return null;
	}

	/**
	 * \brief Returns the Shape object that represents this boundary
	 * 
	 * Returns the Shape object that represents this boundary
	 * 
	 * @return	Shape object that has been constructed to represent this boundary
	 */
	public IsShape getShape() {
		return _myShape;
	}
	
	/**
	 * \brief Return the name of the side of the domain which this boundary is on
	 * 
	 * Return the name of the side of the domain which this boundary is on
	 * 
	 * @return	String containing the name of the side of the domain (e.g. x0z, xNz, etc)
	 */
	public String getSide() {
		return _mySide;
	}
	
	/**
	 * \brief Returns the distance from a point to the boundary
	 * 
	 * Returns the distance from a point to the boundary
	 * 
	 * @param cc	The continuous vector of points to calculate how far the point is from the boundary
	 * @return	Double value stating the distance fromt the point to the boundary
	 */
	public double getDistance(ContinuousVector cc) {
		return _myShape.getDistance(cc);
	}

	/**
	 * \brief For a specified solute, returns the level of that solute in the bulk
	 * 
	 * For a specified solute, returns the level of that solute in the bulk
	 * 
	 * @param soluteIndex	Index of the solute in the simulation dictionary
	 * @return	Value of solute in the connected bulk
	 */
	public double getBulkValue(int soluteIndex) {
		return 0;
	}
	

}
