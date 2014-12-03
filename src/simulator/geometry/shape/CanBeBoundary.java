package simulator.geometry.shape;

import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;

public interface CanBeBoundary extends IsShape
{
	/**
     * Computes orthogonal distance and if this distance is lower than the
     * resolution and if the point is outside, then the point tested is
     * declared to be on the boundary of the domain.
     * 
	 * @param cC	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @param res	Resolution of the domain that this shape is associated
	 * with.
	 * @return	Boolean noting whether this coordinate is on the boundary of
	 * the domain.
	 */
	public Boolean isOnBoundary(ContinuousVector cV, Double res);
	

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
     * \brief Initialisation to create the features of and go along the
     * boundary.
     * 
     * @param aSG	The grid to which this boundary is a part.
     */
	public void readyToFollowBoundary(SpatialGrid aSG);
	
	/**
     * \brief Find the next valid point.
     *  
     *  @param dcIn	Discrete vector within the shape.
     *  @param dcOut	Discrete vector outside the shape.
     *  @param aSG	Spatial grid in which the boundary is associated.
     *  @return Whether a valid point was found.
     */
	public Boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut,
															SpatialGrid aSG);
	
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
	 * \brief Gets the distance from the opposite side (aShape).
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public Double getDistance(CanBeBoundary aBoundary);
	

	
	/**
	 * \brief Return vector normal to the plane.
	 * 
	 * @return	Discrete vector normal to the plane.
	 */
	public DiscreteVector getNormalDC();
}
