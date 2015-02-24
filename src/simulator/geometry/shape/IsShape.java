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

import java.io.Serializable;
import java.util.LinkedList;

import simulator.SpatialGrid;
import simulator.geometry.Domain;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Vertex;
import utils.XMLParser;

/**
 * \brief Interface class used by boundary condition classes. 
 * 
 * Defines the methods used to monitor the boundaries of the computation
 * domain.
 */
public abstract class IsShape implements Serializable
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Indices of the local coordinates used by the Voronoi diagram generator.
	 */
	protected static int _voronoiPrimary, _voronoiSecondary, _voronoiIgnore;
	
	/**
	 * 
	 */
	protected static Double _minPrimary, _maxPrimary;
	
	/**
	 * \brief Reads the coordinates that specify a boundary from the protocol
	 * file, creating a shape.
	 * 
	 * @param shapeRoot	XML elements from the protocol file that contain
	 * coordinates specifying the edge of a boundary.
	 * @param aDomain	The computation domain that this boundary is
	 * associated with.
	 */
	public abstract void readShape(XMLParser shapeRoot, Domain aDomain);
	
	/**
	 * \brief Test if a given set of coordinates is outside this shape.
	 * 
	 * @param cV	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @return	Boolean noting whether this coordinate is inside or outside
	 * this shape.
	 */
	public abstract Boolean isOutside(ContinuousVector cV);
	
	/**
	 * 
	 * @param coord
	 * @return
	 */
	public abstract Boolean isOutside(DiscreteVector coord);
	
	public abstract DiscreteVector getRelativePosition(DiscreteVector coord);
	
	public abstract DiscreteVector getAbsolutePosition(DiscreteVector coord);
	
	public abstract ContinuousVector getRelativePosition(ContinuousVector point);
	
	public abstract ContinuousVector getAbsolutePosition(ContinuousVector point);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape.
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @param ccOut	Corrected coordinates
	 */
	public abstract void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut);
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape,
	 * returning these coordinates.
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @return Corrected coordinates
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		ContinuousVector out = new ContinuousVector();
		orthoProj(ccIn, out);
		return out;
	}
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape.
	 * 
	 * @param dcIn	Coordinates to be corrected
	 * @param dcOut	Corrected coordinates
	 */
	public abstract void orthoProj(DiscreteVector dcIn, DiscreteVector dcOut);
	
	/**
	 * 
	 * @param coord
	 * @return
	 */
	public DiscreteVector getOrthoProj(DiscreteVector coord)
	{
		DiscreteVector out = new DiscreteVector();
		orthoProj(coord, out);
		return out;
	}
	
	/**
	 * \brief Gets the distance from a point on the other side
	 * (ContinuousVector).
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public Double getDistance(ContinuousVector point)
	{
		ContinuousVector diff = getOrthoProj(point);
		diff.subtract(point);
		return diff.norm();
	}
	
	/**
     * \brief Calculates the coordinates of the intersection(s) between a line
     * (point and vector) and the shape.
     * 
     * Intersections must follow the order in which they are encountered,
     * travelling from position along vector.
     * Returns null if none exists. TODO return empty list?
     * 
     * Examples:
     * - Planar may have 0 or 1 intersection
     * - Cylindrical or Hemispherical may have 0, 1, or 2 intersections
     * 
     * @param position	Position used to calculate the line.
     * @param vector	Vector of coordinate positions used to calculate the
     * line.
     * @return Coordinates of the intersection(s) between a line and the shape.
     */
	public abstract LinkedList<ContinuousVector> getIntersections(
						ContinuousVector position, ContinuousVector vector);
	
	
	/**
	 * \brief Takes a vector and returns that vector pointing towards the
	 * inside of the shape.
	 * 
	 * @param cc	Vector outside the shape.
	 * @return ContinuousVector that is pointing towards the inside of the
	 * shape.
	 */
	public abstract ContinuousVector getNormalInside(ContinuousVector cc);
		
	/**
	 * \brief Return vector normal to the plane.
	 * 
	 * @return	Discrete vector normal to the plane.
	 */
	public abstract DiscreteVector getNormalDiscrete();
	
	
	
	/* ------------------------- Voronoi diagram -------------------------- */
	
	public Double getMinPrimary()
	{
		return _minPrimary;
	}
	
	public Double getMaxPrimary()
	{
		return _maxPrimary;
	}
	
	public int getPrimary()
	{
		return _voronoiPrimary;
	}
	
	public int getSecondary()
	{
		return _voronoiSecondary;
	}
	
	/**
	 * \brief gets the 
	 * 
	 * @param point
	 * @return
	 */
	public Double getPrimary(ContinuousVector point)
	{
		Double[] p = convertToLocal(point);
		//System.out.println("p ("+p[0]+", "+p[1]+", "+p[2]+") take "+_voronoiPrimary);
		return p[_voronoiPrimary];
	}
	
	/**
	 * \brief gets the 
	 * 
	 * @param point
	 * @return
	 */
	public Double getSecondary(ContinuousVector point)
	{
		Double[] p = convertToLocal(point);
		return p[_voronoiSecondary];
	}
	
	/**
	 * 
	 * @param site1
	 * @param site2
	 * @return
	 */
	public abstract Edge bisect(Site site1, Site site2);
	
	/**
	 * 
	 * @param point
	 * @return
	 */
	public abstract Double[] convertToLocal(ContinuousVector point);
	
	/**
	 * 
	 * @param local
	 * @return
	 */
	public abstract ContinuousVector convertToVector(Double[] local);
	
	/**
	 * \brief Return the (shortest) distance, over this shape, between two 
	 * points on this shape.
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	public abstract Double distance(ContinuousVector point1, ContinuousVector point2);
	
	/**
	 * 
	 * @param he1
	 * @param he2
	 * @return
	 */
	public abstract Vertex intersect(HalfEdge he1, HalfEdge he2);
	
	
	/* ---------------------- SpatialGrid boundaries ---------------------- */
	
	/**
	 * \brief Gets the distance from the opposite side (aShape).
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public abstract Double getDistance(IsShape aBoundary);
	
	/**
     * \brief Initialisation to create the features of and go along the
     * boundary.
     * 
     * @param aSG	The grid to which this boundary is a part.
     */
	public abstract void readyToFollowBoundary(SpatialGrid aSG);
	
	/**
     * \brief Find the next valid point.
     *  
     *  @param dcIn	Discrete vector within the shape.
     *  @param dcOut	Discrete vector outside the shape.
     *  @param aSG	Spatial grid in which the boundary is associated.
     *  @return Whether a valid point was found.
     */
	public abstract Boolean followBoundary(DiscreteVector dcIn,
									DiscreteVector dcOut, SpatialGrid aSG);
}
