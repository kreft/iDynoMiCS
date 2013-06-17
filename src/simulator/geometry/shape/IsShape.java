/**
 * \package simulator.geometry.shape
 * \brief Package of utilities that assist with managing the boundary conditions in iDynoMiCS
 * 
 * Package of utilities that assist with managing the boundary conditions in iDynoMiCS. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry.shape;

import simulator.geometry.Domain;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.SpatialGrid;
import utils.XMLParser;

/**
 * \brief Interface class used by boundary condition classes. Defines the methods used to monitor the boundaries of the computation domain
 * 
 * Interface class used by boundary condition classes. Defines the methods used to monitor the boundaries of the computation domain
 *
 */
public interface IsShape {

	/**
	 * \brief Reads the coordinates that specify a boundary from the protocol file, creating a shape
	 * 
	 * Reads the coordinates that specify a boundary from the protocol file, creating a shape
	 * 
	 * @param shapeRoot	XML elements from the protocol file that contain coordinates specifying the edge of a boundary
	 * @param aDomain	The computation domain that this boundary is associated with
	 */
	public void readShape(XMLParser shapeRoot,Domain aDomain);
	
	/**
	 * \brief Test if a given set of coordinates is outside this shape
	 * 
	 * Test if a given set of coordinates is outside this shape
	 * 
	 * @param cc	ContinuousVector containing the coordinates of a point to test
	 * @return	Boolean noting whether this coordinate is inside or outside this shape
	 */
	public Boolean isOutside(ContinuousVector cc);
	
	/**
     * Computes orthogonal distance and if this distance is lower than the resolution and if the point is outside, 
     * then the point tested is declared to be on the boundary of the domain
     * 
	 * @param cC	ContinuousVector containing the coordinates of a point to test
	 * @param res	Resolution of the domain that this shape is associated with
	 * @return	Boolean noting whether this coordinate is on the boundary of the domain
	 */
	public Boolean isOnBoundary(ContinuousVector cC, double res);	
	
	/**
     * \brief Calculates the coordinates of the interaction between a line (point a vector) and the plane
     * 
     * Calculates the coordinates of the interaction between a line (point a vector) and the plane. Returns null if none exists
     * 
     * @param position	Position used to calculate the line
     * @param vector	Vector of coordinate positions used to calculate the line
     * @return : coordinates of the intersection between a line and the plane
     */
	public ContinuousVector intersection(ContinuousVector position,ContinuousVector vector);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape 
	 * 
	 * Correct coordinates of a point that has gone outside this shape
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @param ccOut	Corrected coordinates
	 * 
	 */
	public void orthoProj(ContinuousVector ccIn,ContinuousVector ccOut);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape, returning these coordinates
	 * 
	 * Correct coordinates of a point that has gone outside this shape, returning these coordinates
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @return Corrected coordinates
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn);
	
	/**
	 * \brief Used in cyclic boundaries - gets the distance from the opposite side (aShape)
	 * 
	 * Used in cyclic boundaries - gets the distance from the opposite side (aShape)
	 * 
	 * @return Double stating distance to that shape
	 */
	public double getDistance(IsShape aShape);
	
	/**
	 * \brief Used in cyclic boundaries - gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * Used in cyclic boundaries - gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * @return Double stating distance to that shape
	 */
	public double getDistance(ContinuousVector cc);
	
	/**
	 * \brief Takes a vector and returns that vector pointing towards the inside of the shape
	 * 
	 * Takes a vector and returns that vector pointing towards the inside of the shape
	 * 
	 * @param cc	Vector outside the shape
	 * @return ContinuousVector that is pointing towards the inside of the shape
	 * 
	 */
	public ContinuousVector getNormalInside(ContinuousVector cc);
	
	/**
     * \brief Initialisation to create the features of and go along the boundary
     * 
     * Initialisation to create the features of and go along the boundary
     * 
     * @param aSG	The grid to which this boundary is a part
     */
	public void readyToFollowBoundary(SpatialGrid aSG);
	
	/**
     * \brief Find the next valid point
     * 
     *  Find the next valid point 
     *  
     *  @param dcIn	Discrete vector within the shape
     *  @param dcOut	Discrete vector outside the shape
     *  @param aSG	Spatial grid in which the boundary is associated
     *  @return Whether a valid point was found
     *  
     */
	public boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut, SpatialGrid aSG);
	
	/**
	 * \brief Return vector normal to the plane
	 * 
	 * Return vector normal to the plane
	 * 
	 * @return	Discrete vector normal to the plane
	 */
	public DiscreteVector getNormalDC();

}
