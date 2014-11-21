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

import simulator.geometry.*;
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
public class Planar implements IsShape, Serializable 
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
	 * One of two discrete vectors parallel to the plane, i.e. orthogonal to
	 * _vectorDCOut.
	 */
	private DiscreteVector u;
	
	/**
	 * One of two discrete vectors parallel to the plane, i.e. orthogonal to
	 * _vectorDCOut.
	 */
	private DiscreteVector v;
	
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
	public void readShape(XMLParser shapeRoot, Domain aDomain) 
	{
		// Build the variables describing the plane.
		_dPointOnPlane = shapeRoot.getParamIJK("pointIn");
		_dVectorOut = shapeRoot.getParamIJK("vectorOut");
		
		// Translate them into continuous coordinates.
		Double res = aDomain.getGrid().getResolution();
		_cPointOnPlane = new ContinuousVector();
		_cPointOnPlane.x = (_dPointOnPlane.i+(1-_dVectorOut.i)/2)*res;
		_cPointOnPlane.y = (_dPointOnPlane.j+(1-_dVectorOut.j)/2)*res;
		_cPointOnPlane.z = (_dPointOnPlane.k+(1-_dVectorOut.k)/2)*res;
		
		// TODO investigate why the resolution is irrelevant here.
		_cVectorOut = new ContinuousVector();
		_cVectorOut.x = (double) _dVectorOut.i;
		_cVectorOut.y = (double) _dVectorOut.j;
		_cVectorOut.z = (double) _dVectorOut.k;
		
		/* Find two vectors orthogonal to _vectorDCOut, i.e. parallel to the
		 * plane.
		 */
		u = new DiscreteVector();
		v = new DiscreteVector();
		_dVectorOut.orthoVector(u, v);
	}
	
	/**
	 * \brief Test if a given set of coordinates is outside this shape.
	 * 
	 * @param point	ContinuousVector containing the coordinates of a point to
	 * test.
	 * @return	Boolean noting whether this coordinate is inside or outside
	 * this shape.
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		//tempVar.x = -_pointIn.x + point.x;
		//tempVar.y = -_pointIn.y + point.y;
		//tempVar.z = -_pointIn.z + point.z;
		
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
	public Boolean isOnBoundary(ContinuousVector point, Double res)
	{
		return isOutside(point) && (point.distance(getOrthoProj(point))<=res);
	}
	
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
											ContinuousVector vector)
	{
		// If the line is parallel to his plane, return null.
		if ( _cVectorOut.prodScalar(vector).equals(0.0) )
			return null;
		
		// Determine the constant term for the equation of the plane
		Double d = -_cVectorOut.prodScalar(_cPointOnPlane);
				
		
		Double k = - (d + _cVectorOut.prodScalar(position))/
										_cVectorOut.prodScalar(vector);
		
		ContinuousVector out = new ContinuousVector();
		out.x = position.x + k*vector.x;
		out.y = position.y + k*vector.y;
		out.z = position.z + k*vector.z;
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
	public ContinuousVector getNormalInside(ContinuousVector cc)
	{
		return new ContinuousVector(-_cVectorOut.x, -_cVectorOut.y, -_cVectorOut.z);
	}
	
	/**
	 * \brief Correct coordinates of a point that has gone outside this shape. 
	 * 
	 * @param ccIn	Coordinates to be corrected.
	 * @param ccOut	Corrected coordinates.
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		//Double a, b, c, d, k;
		//a = _vectorOut.x;
		//b = _vectorOut.y;
		//c = _vectorOut.z;
		//d = -(_pointIn.x*a + _pointIn.y*b + _pointIn.z*c);
		
		// this next line wasn't calculating the projection coefficient correctly
		// and would case trouble when the plane was away from the substratum
		//k = (d+_vectorOut.prodScalar(ccIn))/(a+b+c);
		// this does it right
		//k = -(d+_vectorOut.prodScalar(ccIn))/_vectorOut.prodScalar(_vectorOut);
		
		// TODO most of these could be calculated just once for the shape
		Double k = _cVectorOut.prodScalar(_cPointOnPlane);
		k -= _cVectorOut.prodScalar(ccIn);
		k /= _cVectorOut.prodScalar(_cVectorOut);
		
		ccOut.x = ccIn.x + k*_cVectorOut.x;
		ccOut.y = ccIn.y + k*_cVectorOut.y;
		ccOut.z = ccIn.z + k*_cVectorOut.z;
	}

	/**
	 * \brief Correct coordinates of a point that has gone outside this shape,
	 * returning these coordinates.
	 * 
	 * @param ccIn	Coordinates to be corrected.
	 * @return Corrected coordinates.
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		ContinuousVector ccOut = new ContinuousVector();
		orthoProj(ccIn, ccOut);
		return ccOut;
	}

	/**
	 * \brief Gets the distance from the opposite side (aShape).
	 * 
	 * Used by cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape.
	 */
	public Double getDistance(IsShape aShape)
	{
		ContinuousVector ccOut = aShape.intersection(_cPointOnPlane, _cVectorOut);
		return _cPointOnPlane.distance(ccOut);
	}
	
	
	/**
	 * \brief Gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * Used in cyclic boundaries.
	 * 
	 * @return Double stating distance to that shape
	 */
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
	public void readyToFollowBoundary(SpatialGrid aSG) 
	{
		Double res = aSG.getResolution();
		
		// Change floor to ceiling?
		origin.i = ((int)Math.floor(_cPointOnPlane.x/res))-u.i;
		origin.j = ((int)Math.floor(_cPointOnPlane.y/res))-u.j;
		origin.k = ((int)Math.floor(_cPointOnPlane.z/res))-u.k;
				
		if( _dVectorOut.i > 0 )
			origin.i--;
		if( _dVectorOut.j > 0 )
			origin.j--;
		if( _dVectorOut.k > 0 )
			origin.k--;
		
		indexU = 0;
		indexV = 0;
		
		if ( u.i == 0 )
			uMax = 0;
		else
			uMax = aSG.getGridSizeI()/u.i;
		if ( u.j != 0 ) 
			uMax = Math.max(uMax, aSG.getGridSizeJ()/u.j);
		if ( u.k != 0 )
			uMax = Math.max(uMax, aSG.getGridSizeK()/u.k);

		if ( v.i == 0 )
			vMax = 0;
		else
			vMax = aSG.getGridSizeI()/v.i;
		if ( v.j != 0 )
			vMax = Math.max(vMax, aSG.getGridSizeJ()/v.j);
		if ( v.k != 0 )
			vMax = Math.max(vMax, aSG.getGridSizeK()/v.k);
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
		move.i = indexU*u.i + indexV*v.i;
		move.j = indexU*u.j + indexV*v.j;
		move.k = indexU*u.k + indexV*v.k;
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
	public DiscreteVector getNormalDC()
	{
		return _dVectorOut;
	}
}