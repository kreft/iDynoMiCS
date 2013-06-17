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

import org.jdom.Element;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.geometry.*;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;

/**
 * \brief Interface for the boundary conditions on the system's margins. To be used in solvers and agent move
 * 
 * Interface for the boundary conditions on the system's margins. To be used in solvers and agent move
 * 
 * @author Andreas Dotsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

public interface IsBC {

	/* Functions related to the geometry of the boundary __________________ */

	/**
     * \brief Used during the initialisation, load the class describing the shape of the boundary defined in the parent class
     * 
     * Used during the initialisation, load the class describing the shape of the boundary defined in the parent class
     * 
     * @param geometryRoot	Usually an XML set of elements that describe the boundary to be created
     * @param aDomain	The computational domain which this boundary is associated with
     */
	public void readGeometry(Element geometryRoot, Domain aDomain);

	/**
	 * \brief Returns the distance from a point to the boundary
	 * 
	 * Returns the distance from a point to the boundary
	 * 
	 * @param cc	The continuous vector of points to calculate how far the point is from the boundary
	 * @return	Double value stating the distance fromt the point to the boundary
	 */
	public double getDistance(ContinuousVector cc);

	/**
     * \brief Calculate the intersection between the boundary and a line (defined by a position and a vector)
     * 
     * Calculate the intersection between the boundary and a line (defined by a position and a vector)
     * 
     * @param position	A continuous vector stating the point to be used in the calculation
     * @param vector	A continuous vector stating the line to be used in the calculation	
     * @return ContinuousVector stating the point of intersection between the boundary and a line
     */
	public ContinuousVector getIntersection(ContinuousVector position, ContinuousVector vector);

	/**
     * \brief Calculate the orthogonal projection of a location on the boundary
     * 
     * Calculate the orthogonal projection of a location on the boundary
     * 
     * @param cc	A continuous vector stating the point to be used in the calculation
     * @return ContinuousVector stating the point on the boundary after the orthogonal projection 
     */
	public ContinuousVector getOrthoProj(ContinuousVector cc);

	/**
	 * \brief Determines if a point is outside the boundary
	 * 
	 * Determines if a point is outside the boundary
	 * 
	 * @param cc	ContinuousVector to check
	 * @return	Boolean value noting whether this coordinate is outside the boundary (true) or not (false)
	 */
	public boolean isOutside(ContinuousVector cc);
	
	/**
	 * \brief Determines if a discrete vector location is outside the boundary
	 * 
	 * Determines if a discrete vector location is outside the boundary
	 * 
	 * @param dC	DiscreteVector to check
	 * @param aSpatialGrid	Spatial grid to check
	 * @return	Boolean value noting whether this coordinate is outside the boundary (true) or not (false)
	 */
	public boolean isOutside(DiscreteVector dC, SpatialGrid aSpatialGrid);
	
	/**
	 * \brief Determine if this boundary is active for solute
	 * 
	 * Determine if this boundary is active for solute
	 * 
	 * @return	Boolean noting whether this boundary is active for solute (true) or not (false)
	 */
	public boolean isActive();
	
	/**
	 * \brief Determine if this boundary is cyclic
	 * 
	 * Determine if this boundary is cyclic
	 * 
	 * @return	Boolean noting whether this boundary is cyclic (true) or not (false)
	 */
	public boolean isCyclic();

	/**
	 * \brief Determine if this boundary is attached to a bulk
	 * 
	 * Determine if this boundary  is attached to a bulk
	 * 
	 * @return	Boolean noting whether this boundary  is attached to a bulk (true) or not (false)
	 */
	public boolean hasBulk();

	
	/**
	 * \brief For a specified solute, returns the level of that solute in the bulk
	 * 
	 * For a specified solute, returns the level of that solute in the bulk
	 * 
	 * @param soluteIndex	Index of the solute in the simulation dictionary
	 * @return	Value of solute in the connected bulk
	 */
	public double getBulkValue(int soluteIndex);
	
	/**
	 * \brief Return the bulk that is connected to this boundary
	 * 
	 * Return the bulk that is connected to this boundary
	 * 
	 * @return Bulk object that is connected to this boundary
	 */
	public Bulk getBulk();
	
	/**
	 * \brief Determine whether this boundary is the supporting structure (substratum)
	 * 
	 * Determine whether this boundary is the supporting structure (substratum)
	 * 
	 * @return Boolean noting whether this boundary is the supporting structure (true) or not (false)
	 */
	public boolean isSupport();

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file. In this case also links the connected bulk to this boundary
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file. In this case also links the connected bulk to this boundary
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBoundCondMarkUp	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, Element aBoundCondMarkUp);

	
	/**
	 * \brief Solver for then boundary condition. Initialises the course along the shape of the boundary during multigrid computation 
	 * 
	 * Solver for the variable concentration boundary condition. Initialises the course along the shape of the boundary during multigrid computation
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid);
	
	/**
	 * \brief Computes and applies gas diffusivity across the gas membrane boundary
	 * 
	 * Computes and applies gas diffusivity across the gas membrane boundary
	 * 
	 * @param relDiff
	 * @param aSolutegrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshDiffBoundary(SoluteGrid relDiff,SoluteGrid aSolutegrid);

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
	public ContinuousVector lookAt(ContinuousVector cc);
	
	/**
     * \brief Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * @param aGroup	LocatedGroup object which has been detected to be outside the boundary
     */
	public void setBoundary(LocatedGroup aGroup);
	
	/**
	 * \brief Modify the movement vector dependent on the boundary condition. In some the agent is killed, in others the agent location is adjusted
	 * 
	 * Modify the movement vector dependent on the boundary condition. In some the agent is killed, in others the agent location is adjusted
	 * 
	 * @param anAgent	The LocatedAgent that has crossed the boundary
	 * @param newLoc	Vector of where this agent was going to be placed
	 */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc);

}