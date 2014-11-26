/**
 * \package simulator.geometry.shape
 * \brief Package of utilities that assist with managing the boundary
 * conditions in iDynoMiCS.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry.shape;

import simulator.geometry.Domain;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.SpatialGrid;
import utils.XMLParser;

/**
 * \brief Interface class used by boundary condition classes. 
 * 
 * Defines the methods used to monitor the boundaries of the computation
 * domain.
 */
public interface IsShape
{
	/**
	 * \brief Reads the coordinates that specify a boundary from the protocol
	 * file, creating a shape.
	 * 
	 * @param shapeRoot	XML elements from the protocol file that contain
	 * coordinates specifying the edge of a boundary.
	 * @param aDomain	The computation domain that this boundary is
	 * associated with.
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain);
	
	/**
	 * \brief Test if a given set of coordinates is outside this shape.
	 * 
	 * @param cV	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @return	Boolean noting whether this coordinate is inside or outside
	 * this shape.
	 */
	public Boolean isOutside(ContinuousVector cV);
	
	/**
     * \brief Calculates the coordinates of the interaction between a line
     * (point a vector) and the plane.
     * 
     * Returns null if none exists.
     * 
     * @param position	Position used to calculate the line.
     * @param vector	Vector of coordinate positions used to calculate the
     * line.
     * @return Coordinates of the intersection between a line and the plane.
     */
	public ContinuousVector intersection(ContinuousVector position,
													ContinuousVector vector);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape.
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @param ccOut	Corrected coordinates
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape,
	 * returning these coordinates.
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @return Corrected coordinates
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn);
	
	/**
	 * \brief Gets the distance from the opposite side (aShape).
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public Double getDistance(IsShape aShape);
	
	/**
	 * \brief Gets the distance from a point on the other side
	 * (ContinuousVector).
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public Double getDistance(ContinuousVector cc);
	
	/**
	 * \brief Takes a vector and returns that vector pointing towards the
	 * inside of the shape.
	 * 
	 * @param cc	Vector outside the shape.
	 * @return ContinuousVector that is pointing towards the inside of the
	 * shape.
	 */
	public ContinuousVector getNormalInside(ContinuousVector cc);
	
	/**
	 * \brief Return vector normal to the plane.
	 * 
	 * @return	Discrete vector normal to the plane.
	 */
	public DiscreteVector getNormalDC();

}
