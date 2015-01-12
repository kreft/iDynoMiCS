/**
 * \package simulator.geometry.shape
 * \brief Package of utilities that assist with managing the boundary
 * conditions in iDynoMiCS.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.geometry.shape;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import simulator.geometry.*;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Vertex;
import simulator.SpatialGrid;
import utils.XMLParser;

/**
 * \brief Create a planar-shaped boundary.
 * 
 * Each boundaryCondition also includes a shape mark-up to define the shape of
 * the boundary. In this release only the 'Planar' class is available, and it
 * requires specification of a point on the boundary and a vector pointing
 * outside the domain. These shape parameters must be given in index 
 * coordinates.
 */
public class Planar implements IsShape, CanBeBoundary, CanPointProcess, Serializable 
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * A point on the plane.
	 */ 
	private DiscreteVector _dPointOnPlane;
	
	/**
	 * A vector normal to the plane and going outside the domain.
	 */
	private DiscreteVector _dVectorOut;
	
	/**
	 * A point on the plane (ContinuousVector).
	 */ 
	private ContinuousVector _cPointOnPlane;
	
	/**
	 * A vector normal to the plan and going outside the domain (ContinuousVector)
	 */
	private ContinuousVector _cVectorOut;
	
	/**
	 * Dot product of the vector out with the point on the plane. Used by
	 * intersection() but only needs to be calculated once. 
	 */
	private Double _vOutDotPPlane;
	
	/**
	 * One of two discrete vectors parallel to the plane, i.e. orthogonal to
	 * _vectorDCOut.
	 */
	private DiscreteVector _dOrthogU;
	
	private ContinuousVector _cOrthogU;
	
	/**
	 * One of two discrete vectors parallel to the plane, i.e. orthogonal to
	 * _vectorDCOut.
	 */
	private DiscreteVector _dOrthogV;
	
	private ContinuousVector _cOrthogV;
	
	/**
	 * Range of discrete coordinates met on this shape
	 */
	private int uMax;
	
	/**
	 * Range of discrete coordinates met on this shape
	 */
	private int vMax;
	
	/**
	 * Index used to check whether a point is within the shape.
	 */
	private static int indexU;
	
	/**
	 * Index used to check whether a point is within the shape.
	 */
	private static int indexV;
	
	/**
	 * Stores the move while being calculated.
	 */
	private static DiscreteVector move = new DiscreteVector();
	
	/**
	 * Temporary store of the origin of a point prior to move.
	 */
	private static DiscreteVector origin = new DiscreteVector();
	
	/**
	 * \brief Reads the coordinates that specify a boundary from the protocol
	 * file, creating a shape.
	 * 
	 * @param shapeRoot	XML elements from the protocol file that contain
	 * coordinates specifying the edge of a boundary.
	 * @param aDomain The computation domain that this boundary is associated
	 * with.
	 */
	@Override
	public void readShape(XMLParser shapeRoot, Domain aDomain) 
	{
		// Build the variables describing the plane.
		_dPointOnPlane = new DiscreteVector(shapeRoot.getParamParser("pointIn"));
		_dVectorOut = new DiscreteVector(shapeRoot.getParamParser("vectorOut"));
		
		/* 
		 * Translate them into continuous coordinates.
		 * TODO investigate why the native 
		 * ContinuousVector(DiscreteVector, res) is not used here 
		 */
		Double res = aDomain.getGrid().getResolution();
		_cPointOnPlane = new ContinuousVector();
		_cPointOnPlane.x = (_dPointOnPlane.i+(1-_dVectorOut.i)/2)*res;
		_cPointOnPlane.y = (_dPointOnPlane.j+(1-_dVectorOut.j)/2)*res;
		_cPointOnPlane.z = (_dPointOnPlane.k+(1-_dVectorOut.k)/2)*res;
		
		// TODO investigate why the resolution is irrelevant here.
		// Try _cVectorOut = new ContinuousVector(_dVectorOut, res); ?
		_cVectorOut = new ContinuousVector(_dVectorOut, res);
		_cVectorOut.x = (double) _dVectorOut.i;
		_cVectorOut.y = (double) _dVectorOut.j;
		_cVectorOut.z = (double) _dVectorOut.k;
		
		_vOutDotPPlane = _cVectorOut.prodScalar(_cPointOnPlane);
		
		/* Find two vectors orthogonal to _vectorDCOut, i.e. parallel to the
		 * plane.
		 */
		_dOrthogU = new DiscreteVector();
		_dOrthogV = new DiscreteVector();
		_dVectorOut.orthoVector(_dOrthogU, _dOrthogV);
		_cOrthogU = new ContinuousVector(_dOrthogU, res);
		_cOrthogU.normalizeVector();
		_cOrthogV = new ContinuousVector(_dOrthogV, res);
		_cOrthogV.normalizeVector();
	}
	
	/**
	 * \brief Test if a given set of coordinates is outside this shape.
	 * 
	 * @param point	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @return	Boolean noting whether this coordinate is inside or outside
	 * this shape.
	 */
	@Override
	public Boolean isOutside(ContinuousVector point)
	{
		ContinuousVector temp = new ContinuousVector(point);
		temp.subtract(_cPointOnPlane);
		return ( _cVectorOut.cosAngle(temp) > 0 );
	}
	
	/**
     * Computes orthogonal distance and if this distance is lower than the
     * resolution and if the point is outside, then the point tested is
     * declared to be on the boundary of the domain.
     * 
	 * @param point	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @param res	Resolution of the domain that this shape is associated
	 * with.
	 * @return	Boolean noting whether this coordinate is on the boundary of
	 * the domain.
	 */
	@Override
	public Boolean isOnBoundary(ContinuousVector point, Double res)
	{
		return isOutside(point) && (point.distance(getOrthoProj(point))<=res);
	}
	
	/**
     * \brief Calculates the coordinates of the interaction between a line 
     * (point and vector) and the plane.
     * 
     * Returns null if none exists.
     * 
     * @param position	Position used to calculate the line.
     * @param vector	Vector of coordinate positions used to calculate the
     * line.
     * @return Coordinates of the intersection between a line and the plane.
     */
	@Override
	public ContinuousVector intersection(ContinuousVector position, 
											ContinuousVector vector)
	{
		// If the line is parallel to his plane, return null.
		Double c = _cVectorOut.prodScalar(vector);
		if ( c.equals(0.0) )
			return null;
		ContinuousVector out = new ContinuousVector(vector);
		/* Find the (relative) length along vector we must travel to hit the
		 * plane.
		 */
		out.times((_vOutDotPPlane - _cVectorOut.prodScalar(position))/c);
		out.add(position);
		return out;
	}
	
	/**
	 * \brief Takes a vector and returns that vector pointing towards the
	 * inside of the shape.
	 * 
	 * @param cc	Vector outside the shape.
	 * @return ContinuousVector that is pointing towards the inside of the
	 * shape.
	 * 
	 */
	@Override
	public ContinuousVector getNormalInside(ContinuousVector cc)
	{
		ContinuousVector out = new ContinuousVector(_cVectorOut);
		out.turnAround();
		return out;
	}
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape. 
	 * 
	 * @param ccIn	Coordinates to be corrected.
	 * @param ccOut	Corrected coordinates.
	 */
	@Override
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		ccOut = intersection(ccIn, _cVectorOut);
	}

	/**
	 * \brief Correct coordinates of a point that has gone outside this shape,
	 * returning these coordinates.
	 * 
	 * @param ccIn	Coordinates to be corrected.
	 * @return Corrected coordinates.
	 * 
	 */
	@Override
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		return intersection(ccIn, _cVectorOut);
	}

	/**
	 * \brief Gets the distance from the opposite side (aShape).
	 * 
	 * Used by cyclic boundaries.
	 * 
	 * TODO Rob Dec2014: This is only guaranteed to be correct if planes are
	 * parallel.
	 * 
	 * @return Double stating distance to that shape.
	 */
	@Override
	public Double getDistance(CanBeBoundary aBoundary)
	{
		ContinuousVector ccOut = aBoundary.intersection(_cPointOnPlane, _cVectorOut);
		return _cPointOnPlane.distance(ccOut);
	}
	
	
	/**
	 * \brief Gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape
	 */
	@Override
	public Double getDistance(ContinuousVector cc)
	{
		ContinuousVector ccOut = intersection(cc, _cVectorOut);
		return ccOut.distance(cc);
	}
	
	/**
     * \brief Initialisation to create the features of and go along the
     * boundary.
     * 
     * @param aSG	The grid to which this boundary is a part.
     */
	@Override
	public void readyToFollowBoundary(SpatialGrid aSG) 
	{
		Double res = aSG.getResolution();
		
		// Change floor to ceiling?
		origin.i = ((int)Math.floor(_cPointOnPlane.x/res))-_dOrthogU.i;
		origin.j = ((int)Math.floor(_cPointOnPlane.y/res))-_dOrthogU.j;
		origin.k = ((int)Math.floor(_cPointOnPlane.z/res))-_dOrthogU.k;
				
		if( _dVectorOut.i > 0 )
			origin.i--;
		if( _dVectorOut.j > 0 )
			origin.j--;
		if( _dVectorOut.k > 0 )
			origin.k--;
		
		indexU = 0;
		indexV = 0;
		
		if ( _dOrthogU.i == 0 )
			uMax = 0;
		else
			uMax = aSG.getGridSizeI()/_dOrthogU.i;
		if ( _dOrthogU.j != 0 ) 
			uMax = Math.max(uMax, aSG.getGridSizeJ()/_dOrthogU.j);
		if ( _dOrthogU.k != 0 )
			uMax = Math.max(uMax, aSG.getGridSizeK()/_dOrthogU.k);

		if ( _dOrthogV.i == 0 )
			vMax = 0;
		else
			vMax = aSG.getGridSizeI()/_dOrthogV.i;
		if ( _dOrthogV.j != 0 )
			vMax = Math.max(vMax, aSG.getGridSizeJ()/_dOrthogV.j);
		if ( _dOrthogV.k != 0 )
			vMax = Math.max(vMax, aSG.getGridSizeK()/_dOrthogV.k);
	}

	/**
     * \brief Find the next valid point.
     * 
     *  @param dcIn	Discrete vector within the shape
     *  @param dcOut	Discrete vector outside the shape
     *  @param aSG	Spatial grid in which the boundary is associated
     *  @return Whether a valid point was found
     *  
     */
	@Override
	public Boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut,
															SpatialGrid aSG) 
	{
		// Find the next valid point
		Boolean vectorValid = false;
		do 
		{
			stepBoundary();
			dcIn.sendSum(origin, move);
			vectorValid = aSG.isValid(dcIn);
			// If a valid point has been found, compute its closest neighbour outside
			if (vectorValid)
				dcOut.sendSum(dcIn, _dVectorOut);
		} while ( !(vectorValid) && indexV<vMax );
		
		return vectorValid; 
	}

	/**
     * \brief Process next location on the boundary.
     */
	public void stepBoundary()
	{
		if ( indexU < uMax )
			indexU++;
		else
		{
			indexU = 0;
			indexV++;
		}
		move.i = indexU*_dOrthogU.i + indexV*_dOrthogV.i;
		move.j = indexU*_dOrthogU.j + indexV*_dOrthogV.j;
		move.k = indexU*_dOrthogU.k + indexV*_dOrthogV.k;
	}

	/**
     * \brief Test that the received point is coplanar with the definition
     * point.
     * 
     * @param point	ContinuousVector of coordinates that should be tested
     * @return	Boolean stating whether the point is cooplanar with the definition point
     */
	public boolean isOnShape(ContinuousVector point)
	{
		Double xDiff = (point.x-_cPointOnPlane.x)/_cVectorOut.x; 
		boolean out = (xDiff == (point.y-_cPointOnPlane.y)/_cVectorOut.y);
		out &= (xDiff == (point.z-_cPointOnPlane.z)/_cVectorOut.z);
		return out;
	}

	
	/**
	 * \brief Return vector normal to the plane.
	 * 
	 * @return	Discrete vector normal to the plane
	 */
	@Override
	public DiscreteVector getNormalDC()
	{
		return _dVectorOut;
	}
	
	/**
	 * Currently assumes both points on the plane.
	 */
	@Override
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		return point1.distance(point2);
	}
	
	/**
	 * TODO check
	 */
	@Override
	public Edge bisect(Site site1, Site site2)
	{
		if ( site1.equals(site2) )
			return null;
		
		Edge out = new Edge();
		
		out.region[0] = site1;
		out.region[1] = site2;
		
		// Already implicitly the case, but stated explicitly for clarity.
		out.endPoint[0] = null;
		out.endPoint[1] = null;
		
		ContinuousVector difference = new ContinuousVector(site2);
		difference.subtract(site1);
		
		// TODO change multiplier to account for site weighting
		Double weight = 0.5;
		
		ContinuousVector midPoint = new ContinuousVector(difference);
		midPoint.times(weight);
		midPoint.add(site1);
		
		/* Find the difference and the midPoint in U, V coordinates. Since
		 * _cOrthogU and _cOrthogV have been normailsed, this is equivalent to
		 * dU = difference.norm() * _cOrthogU.cosAngle(difference);
		 * etc
		 */
		Double dU = _cOrthogU.prodScalar(difference);
		Double dUmid = _cOrthogU.prodScalar(midPoint);
		Double dV = _cOrthogV.prodScalar(difference);
		Double dVmid = _cOrthogV.prodScalar(midPoint);
		
		// c0*U + c1*V = c2 
		if ( Math.abs(dU) > Math.abs(dV) )
		{
			// u + slope*v = midU + slope*midV
			out.coefficient[0] = 1.0;
			out.coefficient[1] = dV/dU;
			// Don't need to worry about dU = 0, as |dU| > |dV|
			out.coefficient[2] = dUmid + (dVmid*dV/dU);
		}
		else
		{
			// u/slope + v = midU/slope + midV
			out.coefficient[0] = dU/dV;
			out.coefficient[1] = 1.0;
			// Don't need to worry about dV = 0, as |dU| < |dV|
			out.coefficient[2] = (dUmid*dU/dV) + dVmid;
		}
		
		return out;
	}
	
	public Vertex intersect(HalfEdge he1, HalfEdge he2)
	{
		Vertex out = new Vertex();
		// TODO
		return out;
	}
	
	/**
	 * 
	 * 
	 *
	 */
	public final int compare(ContinuousVector point1,
												ContinuousVector point2)
	{		
		int out = (int) Math.signum(point1.y - point2.y);
		if ( out == 0 )
			out = (int) Math.signum(point1.x - point2.x);
		return out;

	}
}
