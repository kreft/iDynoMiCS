/**
 * \package simulator.geometry.boundaryConditions
 * \brief Package of boundary conditions that can be used to capture agent
 * behaviour at the boundary of the computation domain.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.geometry.boundaryConditions;

import java.util.LinkedList;

import utils.LogFile;
import utils.XMLParser;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import simulator.geometry.shape.*;

/**
 * \brief Group all methods expected by the interface but common to most of
 * the boundary classes.
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public abstract class AllBC
{
	/* _____________________________ FIELDS _______________________________ */
	/**
	 * The name of the boundary describing its side (xOy,...).
	 */
	protected String	_mySide;
	
	/**
	 * The shape of the boundary.
	 */
	protected IsShape	_myShape;
	
	/**
	 * Boolean noting whether this boundary is the supporting structure
	 * (substratum).
	 */
	protected boolean	_isSupport	= false;
	
	/**
	 * Boolean noting whether this boundary can contain active solute.
	 */
	protected boolean	activeForSolute	= true;
	
	/* _________________ INTERNAL TEMPRARY VARIABLES ______________________ */
	/**
	 * Discrete coordinates of a voxel inside the computation domain but along
	 * the boundary.
	 */ 
	protected static DiscreteVector dcIn	= new DiscreteVector();
	
	/**
	 * Discrete coordinates of the voxel in front of the one outside the boundary
	 */
	protected static DiscreteVector dcOut	= new DiscreteVector();

	/* ________________________ CONSTRUCTION METHODS ________________________ */

	/**
	 * \brief Generic constructor called to dynamically instantiate a child
	 * class object.
	 * 
	 * @param root	Set of XML tags relating to one boundary condition
	 * @param aSim	The current simulation object used to simulate the
	 * conditions specified in this protocol file.
	 * @param aDomain	The computation domain to which this boundary
	 * condition is assigned.
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
	 * \brief Initialises the boundary condition.
	 * 
	 * This method should be overridden by each boundary condition class file.
	 * 
	 * @param aSim	The current simulation object used to simulate the
	 * conditions specified in this protocol file.
	 * @param aDomain	The computation domain to which this boundary
	 * condition is assigned.
	 * @param root	Set of XML tags relating to one boundary condition.
	 */
	public abstract void init(Simulator aSim, Domain aDomain, XMLParser root);
	
	/**
     * \brief Used during the initialisation, load the class describing the
     * shape of the boundary defined in the parent class.
     * 
     * @param geometryRoot	Usually an XML set of elements that describe the
     * boundary to be created.
     * @param aDomain	The computational domain which this boundary is
     * associated with.
     */
	public void readGeometry(XMLParser geometryRoot, Domain aDomain)
	{
		// Set the name of the boundary.
		_mySide = geometryRoot.getName();
		// Set the class to use to define the shape.
		String className = "simulator.geometry.shape.";
		className += geometryRoot.getChildParser("shape").getAttribute("class");
		// Build the instance used to describe the shape.
		try
		{
			_myShape = (IsShape) Class.forName(className).newInstance();
			_myShape.readShape(new 
					XMLParser(geometryRoot.getChildElement("shape")), aDomain);
		}
		catch (Exception e)
		{
			
		}
	}
	
	/**
	 * \brief Determines if a point is outside the boundary.
	 * 
	 * @param position	ContinuousVector to check.
	 * @return	Boolean value noting whether this coordinate is outside the
	 * boundary (true) or not (false).
	 */
	public Boolean isOutside(ContinuousVector position)
	{
		return _myShape.isOutside(position);
	}
	
	/**
	 * \brief Return the name of the side of the domain which this boundary is
	 * on.
	 * 
	 * TODO this is a duplicate of getSide()!
	 * 
	 * @return	String containing the name of the side of the domain
	 * (e.g. x0z, xNz, etc).
	 */
	public String getSideName()
	{
		return this._mySide;
	}
	
	/**
	 * \brief Solver for the variable concentration boundary condition.
	 *  
	 * Initialises the course along the shape of the boundary during multigrid
	 * computation.
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be
	 * refreshed by the solver.
	 * @see ComputationDomain.refreshBoundaries()
	 */
	public abstract void refreshBoundary(SoluteGrid aSoluteGrid);
	
	/**
     * \brief Method used if a boundary modifies the local diffusivity
     * constant. Most of boundaries do not modify it.
     * 
     * @param relDiff	Relative difference grid
     * @param aSolutegrid	Grid of solute information which is to be
     * refreshed by the solver.
     * 
     * @see BoundaryGasMembrane
     */
	public void refreshDiffBoundary(SoluteGrid relDiff, SoluteGrid aSolutegrid)
	{
		
	};

	/* ______________INTERACTION WITH THE PARTICLES _____________________ */

	/**
	 * \brief Method used by another which gets the indexed grid position of a
	 * continuous vector.
	 * 
	 * Some boundary conditions (e.g. BoundaryCyclic) need the input corrected
	 * due to the condition, some don't and just return the input. Maybe we'll
	 * change this at some point as to just return the input looks a bit daft
	 * - but we'll leave it here for the moment.
	 * 
	 * @param position ContinuousVector that gives the current location of an
	 * agent to check on the grid.
	 */
	public ContinuousVector lookAt(ContinuousVector position)
	{
		return position;
	}
	
	/**
     * \brief Change the status of a specified LocatedGroup to note that it
     * has been identified as being outside this boundary.
     * 
     * @param aGroup	LocatedGroup object which has been detected to be
     * outside the boundary.
     */
	public abstract void setBoundary(LocatedGroup aGroup);
	
	/**
     * Modify the movement vector : the new position is the orthogonal 
     * projection on the boundary surface.
     * 
     * @see LocatedAgent.move();
     */
	public abstract void applyBoundary(LocatedAgent anAgent, 
													ContinuousVector newLoc);
	
	public void applyBoundary(DiscreteVector coord)
	{
		if ( _myShape.isOutside(coord) )
			coord = _myShape.getOrthoProj(coord);
	}
	
	/* ___________________ INTERACTION WITH THE DOMAIN _________________________ */

	/**
	 * \brief Determine whether this boundary is the supporting structure
	 * (substratum).
	 * 
	 * @return Boolean noting whether this boundary is the supporting
	 * structure (true) or not (false).
	 */
	public boolean isSupport()
	{
		return _isSupport;
	}

	/**
	 * \brief Determine if this boundary is active for solute.
	 * 
	 * @return	Boolean noting whether this boundary is active for solute
	 * (true) or not (false).
	 */
	public boolean isActive()
	{
		return activeForSolute;
	}
	
	/**
	 * \brief Determines if a discrete vector location is outside the boundary.
	 * 
	 * @param dc	DiscreteVector to check.
	 * @param aSpatialGrid	The grid to check whether a point is outside.
	 * @return	Boolean value noting whether this coordinate is outside the
	 * boundary (true) or not (false).
	 */
	public boolean isOutside(DiscreteVector dc, SpatialGrid aSpatialGrid)
	{
		ContinuousVector temp = new ContinuousVector();
		temp.setToVoxelCenter(dc, aSpatialGrid.getResolution());
		return _myShape.isOutside(temp);
	}

	/* ____________________ TOOLBOX ______________________________ */

	/**
     * \brief Calculate the intersection(s) between the boundary and a line
     * (defined by a position and a vector).
     * 
     * @param position	A continuous vector stating the point to be used in
     * the calculation.
     * @param vector	A continuous vector stating the line to be used in the
     * calculation.	
     * @return List of ContinuousVectors stating the point(s) of intersection
     * between the boundary and a line.
     */
	public LinkedList<ContinuousVector> getIntersections(
						ContinuousVector position, ContinuousVector vector)
	{
		return _myShape.getIntersections(position, vector);
	}
	
	/**
     * \brief Calculate the orthogonal projection of a location on the
     * boundary.
     * 
     * @param position	A continuous vector stating the point to be used in
     * the calculation.
     * @return ContinuousVector stating the point on the boundary after the
     * orthogonal projection. 
     */
	public ContinuousVector getOrthoProj(ContinuousVector position)
	{
		return _myShape.getOrthoProj(position);
	}

	/**
	 * \brief Returns the Shape object that represents this boundary.
	 * 
	 * @return	Shape object that has been constructed to represent this
	 * boundary.
	 */
	public IsShape getShape()
	{
		return _myShape;
	}
	
	/**
	 * \brief Return the name of the side of the domain which this boundary is
	 * on.
	 * 
	 * @return	String containing the name of the side of the domain (e.g.
	 * x0z, xNz, etc).
	 */
	public String getSide()
	{
		return _mySide;
	}
	
	/**
	 * \brief Returns the distance from a point to the boundary.
	 * 
	 * @param position	The continuous vector of points to calculate how far 
	 * the point is from the boundary.
	 * @return	Double value stating the distance from the point to the
	 * boundary.
	 */
	public Double getDistance(ContinuousVector position)
	{
		return _myShape.getDistance(position);
	}
	
	/**
	 * \brief Kills an agent trying to leave this boundary.
	 * 
	 * Called by subclass boundaries.
	 * 
	 * @param anAgent
	 * @param target
	 * @param reason
	 * 
	 * @see applyBoundary(LocatedAgent anAgent, ContinuousVector target)
	 */
	protected void deadlyBoundary(LocatedAgent anAgent,
									ContinuousVector target, String reason)
	{
		/*
		 * Recording reason of death: agent will be moved to agentToKill list
		 * when die() calls registerDeath().
		 */
		anAgent.death = reason;
		
		anAgent.die(false);
		// To label this agent as "shoving solved", set to zero its movement.
		anAgent.getMovement().reset();
		target.set(anAgent.getLocation());
	}
	
	/**
	 * \brief Stops an agent trying to leave this boundary.
	 * 
	 * Called by subclass boundaries.
	 * 
	 * @param anAgent
	 * @param target
	 * 
	 * @see applyBoundary(LocatedAgent anAgent, ContinuousVector target)
	 */
	protected void hardBoundary(LocatedAgent anAgent, ContinuousVector target)
	{
		// Define coordinates of the corrected position.
		_myShape.orthoProj(target, target);
		/*
		 * Build a vector normal to the boundary and starting from the
		 * orthogonal projection.
		 */
		ContinuousVector vectorIn = 
							new ContinuousVector(_myShape.getNormalInside());
		/*
		 * The whole cell has to be inside, so make a translation equal to the
		 * total radius.
		 */
		vectorIn.normalizeVector();
		vectorIn.times(anAgent.getRadius(true));
		// Compute the new position.
		target.add(vectorIn);
		// Compute and update the movement vector leading to this new position.
		anAgent.getMovement().sendDiff(target,anAgent.getLocation());
	}
	
	public Double overBoundary(Double radius, ContinuousVector target) {
		Double v = _myShape.getDistance(target)-radius;
		if (isOutside(target) || (v < 0))
			return v;
		else
			return null;

		
	}
}
