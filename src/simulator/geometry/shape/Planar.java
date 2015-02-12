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

import java.util.LinkedList;

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
public class Planar extends IsShape 
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
	 * Dot product of the vector out with the point on the plane (continuous).
	 * Used by intersection() but only needs to be calculated once. 
	 */
	private Double _cOutDotPPlane;
	
	/**
	 * Dot product of the vector out with the point on the plane (discrete).
	 */
	private int _dOutDotPPlane;
	
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
	 * \brief Constructor method of a planar shape.
	 * 
	 * Used in testing.
	 * 
	 * @param dPointOnPlane Discrete vector giving the position of a point on
	 * the plane.
	 * @param dVectorOut Discrete vector giving the direction of a vector
	 * orthogonal (i.e. at right angles to) the plane.
	 * @param resolution The resolution of the domain grid.
	 */
	public Planar(DiscreteVector dPointOnPlane,
								DiscreteVector dVectorOut, Double resolution)
	{
		_dPointOnPlane = new DiscreteVector(dPointOnPlane);
		_dVectorOut = new DiscreteVector(dVectorOut);
		init(resolution);
	}
	
	public void readShape(XMLParser shapeRoot, Domain aDomain) 
	{
		// Build the variables describing the plane.
		_dPointOnPlane = new DiscreteVector(shapeRoot.getParamParser("pointIn"));
		_dVectorOut = new DiscreteVector(shapeRoot.getParamParser("vectorOut"));
		
		init(aDomain.getResolution());
	}
	
	/**
	 * 
	 * @param resolution
	 */
	private void init(Double resolution)
	{
		/*
		 * TODO investigate why the native 
		 * ContinuousVector(DiscreteVector, res) is not used here.
		 */
		_cPointOnPlane = new ContinuousVector();
		_cPointOnPlane.x = (_dPointOnPlane.i+(1-_dVectorOut.i)/2)*resolution;
		_cPointOnPlane.y = (_dPointOnPlane.j+(1-_dVectorOut.j)/2)*resolution;
		_cPointOnPlane.z = (_dPointOnPlane.k+(1-_dVectorOut.k)/2)*resolution;
		// The vector out just needs to have the correct direction.
		_cVectorOut = new ContinuousVector();
		_cVectorOut.x = new Double(_dVectorOut.i);
		_cVectorOut.y = new Double(_dVectorOut.j);
		_cVectorOut.z = new Double(_dVectorOut.k);
		/*
		 * Dots products of the vector out with the point on the plane. Saves
		 * having to calculate these values repeatedly.
		 */
		_cOutDotPPlane = _cVectorOut.prodScalar(_cPointOnPlane);
		_dOutDotPPlane = _dVectorOut.prodScalar(_dPointOnPlane);
		/* 
		 * Find two vectors orthogonal to each vector out, i.e. parallel to the
		 * plane.
		 */
		_dOrthogU = new DiscreteVector();
		_dOrthogV = new DiscreteVector();
		_dVectorOut.orthoVector(_dOrthogU, _dOrthogV);
		
		_cOrthogU = new ContinuousVector();
		_cOrthogU.set(_dOrthogU);
		_cOrthogU.normalizeVector();
		
		_cOrthogV = new ContinuousVector();
		_cOrthogV.set(_dOrthogV);
		_cOrthogV.normalizeVector();
		
		/*
		System.out.println("discrete norm: "+_dVectorOut.toString());
		System.out.println("u: "+_dOrthogU.toString());
		System.out.println("v: "+_dOrthogV.toString());
		System.out.println("continuous norm: "+_cVectorOut.toString());
		System.out.println("u: "+_cOrthogU.toString());
		System.out.println("v: "+_cOrthogV.toString());
		*/
		
		_voronoiPrimary = 0;
		_voronoiSecondary = 1;
	}
	
	public Boolean isOutside(ContinuousVector point)
	{
		return ( _cVectorOut.cosAngle(getRelativePosition(point)) > 0 );
	}
	
	public Boolean isOutside(DiscreteVector coord)
	{
		return ( _dVectorOut.cosAngle(getRelativePosition(coord)) > 0 );
	}
	
	public LinkedList<ContinuousVector> getIntersections(
						ContinuousVector position, ContinuousVector vector)
	{
		// If the line is parallel to his plane, return null.
		Double c = _cVectorOut.prodScalar(vector);
		if ( c.equals(0.0) )
			return null;
		LinkedList<ContinuousVector> out = new LinkedList<ContinuousVector>();
		ContinuousVector intersection = new ContinuousVector(vector);
		/* 
		 * Find the (relative) length along vector we must travel to hit the
		 * plane.
		 */
		intersection.times((_cOutDotPPlane - _cVectorOut.prodScalar(position))/c);
		intersection.add(position);
		out.add(intersection);
		return out;
	}
	
	/**
	 * 
	 * @param point
	 * @return
	 */
	public ContinuousVector getRelativePosition(ContinuousVector point)
	{
		ContinuousVector pointOnPlaneToPoint = new ContinuousVector();
		pointOnPlaneToPoint.sendDiff(point, _cPointOnPlane);
		return pointOnPlaneToPoint;
	}
	
	public ContinuousVector getAbsolutePosition(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.add(_cPointOnPlane);
		return out;
	}
	
	public DiscreteVector getRelativePosition(DiscreteVector coord)
	{
		DiscreteVector pointOnPlaneToPoint = new DiscreteVector();
		pointOnPlaneToPoint.sendDiff(coord, _dPointOnPlane);
		return pointOnPlaneToPoint;
	}
	
	public DiscreteVector getAbsolutePosition(DiscreteVector coord)
	{
		DiscreteVector out = new DiscreteVector();
		out.sendSum(coord, _dPointOnPlane);
		return out;
	}
	
	public Double[] convertToLocal(ContinuousVector point)
	{
		Double[] out = new Double[3];
		ContinuousVector pointOnPlaneToPoint = getRelativePosition(point);
		out[0] = _cOrthogU.cosAngle(pointOnPlaneToPoint);
		out[1] = _cOrthogV.cosAngle(pointOnPlaneToPoint);
		out[2] = _cVectorOut.cosAngle(pointOnPlaneToPoint);
		return out;
	}
	
	public ContinuousVector convertToVector(Double[] local)
	{
		ContinuousVector out = new ContinuousVector();
		ContinuousVector temp = new ContinuousVector(_cVectorOut);
		temp.times(local[0]);
		out.add(temp);
		temp.set(_cOrthogV);
		temp.times(local[1]);
		out.add(temp);
		temp.set(_cVectorOut);
		temp.times(local[2]);
		out.add(temp);
		out.add(_cPointOnPlane);
		return out;
	}
	
	public ContinuousVector getNormalContinuous()
	{
		return _cVectorOut;
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
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		ccOut = getIntersections(ccIn, _cVectorOut).getFirst();
	}
	
	@Override
	public void orthoProj(DiscreteVector dcIn, DiscreteVector dcOut)
	{
		/*
		 * TODO Rob 23Jan2015: This is not as robust as I would like!
		 * It will work fine if _dVector out is parallel to one axis only
		 * (as is usually the case) but may fail otherwise.
		 */
		dcOut.set(_dVectorOut);
		dcOut.times(_dOutDotPPlane - _dVectorOut.prodScalar(dcIn));
		dcOut.add(dcIn);
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
	public Double getDistance(IsShape aBoundary)
	{
		Double out = Double.MAX_VALUE;
		for ( ContinuousVector ccOut : 
					aBoundary.getIntersections(_cPointOnPlane, _cVectorOut) )
		{
			out = Math.min(out, _cPointOnPlane.distance(ccOut));
		}
		return out;
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
		ContinuousVector ccOut = getIntersections(cc, _cVectorOut).getFirst();
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
	 * \brief Return vector normal to the plane.
	 * 
	 * @return	Discrete vector normal to the plane
	 */
	@Override
	public DiscreteVector getNormalDiscrete()
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
		ContinuousVector direction = _cVectorOut.crossProduct(difference);
		direction.normalizeVector();
		
		ContinuousVector midPoint = new ContinuousVector(difference);
		midPoint.times(0.5);
		midPoint.add(site1);
		
		return lineToEdge(midPoint, direction);
	}
	
	public Vertex intersect(Edge edge1, Edge edge2)
	{
		Vertex out = new Vertex();
		
		return out;
	}
	
	public Vertex intersect(HalfEdge he1, HalfEdge he2)
	{
		return intersect(he1.edge, he2.edge);
	}
	
	public Edge intersect(Planar other) 
	{
		// Check planes are not parallel.
		if ( _dVectorOut.isParallel(other.getNormalDiscrete()) )
			return null;
		// The edge will be parallel to the cross product of the two norms.
		ContinuousVector direction = 
						_cVectorOut.crossProduct(other.getNormalContinuous());
		System.out.println("dir: "+direction.toString());
		// Find a point on the Edge.
		ContinuousVector point = direction.crossProduct(_cVectorOut);
		System.out.println("point0: "+point.toString());
		point = other.getIntersections(_cPointOnPlane, point).getFirst();
		System.out.println("point1: "+point.toString());
		return lineToEdge(point, direction);
	}
	
	/**
	 * TODO Check thoroughly!
	 * 
	 * @param point
	 * @param direction
	 * @return
	 */
	public Edge lineToEdge(ContinuousVector point, ContinuousVector direction)
	{
		Double[] p = convertToLocal(point);
		System.out.println("p: ("+p[0]+","+p[1]+","+p[2]+")");
		Double[] d = convertToLocal(direction);
		System.out.println("d: ("+d[0]+","+d[1]+","+d[2]+")");
		/*
		 * Check the point is on the plane, that the line is parallel to it,
		 * and that the direction vector is non-zero.   
		 */
		if ( p[2] != 0.0 || d[2] != 0.0 || 
				(d[_voronoiPrimary] == 0.0 && d[_voronoiSecondary] == 0.0) )
		{
			return null;
		}
		Edge out = new Edge();
		if ( Math.abs(d[_voronoiPrimary]) > Math.abs(d[_voronoiSecondary]) )
		{
			out.coefficient[0] = 1.0;
			out.coefficient[1] = d[_voronoiSecondary]/d[_voronoiPrimary];
			out.coefficient[2] = p[_voronoiPrimary] + (p[1]*out.coefficient[1]);
		}
		else
		{
			out.coefficient[0] = d[_voronoiPrimary]/d[_voronoiSecondary];
			out.coefficient[1] = 1.0;
			out.coefficient[2] = (p[_voronoiPrimary]*out.coefficient[0]) + p[_voronoiSecondary];
		}
		return out;
	}
	
	/**
	 * \brief Converts a HalfEdge back to a pair of ContinuousVectors.
	 * 
	 * The first vector is gives a point on the line, the second the direction
	 * of the line. The direction vector may not have the same length as the
	 * original direction vector, and may even be pointing the opposite way,
	 * but it will be parallel.
	 * 
	 * @param edge
	 * @return
	 */
	public ContinuousVector[] edgeToLine(Edge edge)
	{
		Double[] point = new Double[3];
		Double[] direction = new Double[3];
		point[2] = 0.0;
		direction[2] = 0.0;
		if ( edge.coefficient[0] == 1.0 )
		{
			/* This corresponds to
			 * Math.abs(d[_voronoiPrimary]) > Math.abs(d[_voronoiSecondary])
			 * in lineToHalfEdge(point, direction)
			 */
			point[_voronoiPrimary] = edge.coefficient[2] - edge.coefficient[1];
			point[_voronoiSecondary] = 1.0;
			direction[_voronoiPrimary] = 1.0;
			direction[_voronoiSecondary] = edge.coefficient[1];
		}
		else
		{
			point[0] = 1.0;
			point[1] = edge.coefficient[2] - edge.coefficient[1];
			direction[0] = edge.coefficient[1];
			direction[1] = 1.0;
		}
		ContinuousVector[] out = new ContinuousVector[2];
		out[0] = convertToVector(point);
		out[1] = convertToVector(direction);
		return out;
	}
	
	public void restrictPlane(LinkedList<Planar> walls)
	{
		_maxStar = -Double.MAX_VALUE;
		_minStar = Double.MAX_VALUE;
		LinkedList<Edge> limits = new LinkedList<Edge>();
		Edge temp;
		Double[] v;
		for ( Planar plane : walls )
		{
			System.out.println("Looking at "+plane.getNormalDiscrete().toString());
			temp = intersect(plane);
			if ( temp != null )
				limits.add(temp);
		}
		for ( Edge limit : limits )
			for ( Edge other : limits )
			{
				if ( limit.equals(other) )
					continue;
				v = convertToLocal(intersect(limit, other));
				_maxStar = Math.max(_maxStar, v[_voronoiSecondary]);
				_minStar = Math.min(_minStar, v[_voronoiSecondary]);
			}
	}
}
