package simulator.geometry.shape;

import java.util.LinkedList;

import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.Domain;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Vertex;
import utils.Complex;
import utils.ExtraMath;
import utils.XMLParser;

public class Cylindrical extends IsShape
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * (0) height, (1) radius, (2) angle
	 */
	public final static int azimuthCoord = 2;
	
	/**
	 * 
	 */
	public final static int heightCoord = 0;
	
	/**
	 * 
	 */
	public final static int radialCoord = 1;
	
	/**
	 * A point on the cylinder axis
	 */
	private DiscreteVector _dPointCenterBase;
	
	/**
	 * 
	 */
	private DiscreteVector _dVectorAlongAxis;
	
	/**
	 * A vector, orthogonal to the vector along the axis, used as a reference
	 * for the angle theta.
	 */
	private DiscreteVector _dVectorRadiusV;
	
	/**
	 * Another vector, orthogonal to the vector along the axis, used for
	 * converting back to Cartesian coordinates.
	 */
	private DiscreteVector _dVectorRadiusW;
	
	private ContinuousVector _cPointCenterBase;
	
	
	private ContinuousVector _cVectorAlongAxis;
	
	/**
	 * A vector, orthogonal to the vector along the axis, used as a reference
	 * for the angle theta.
	 */
	private ContinuousVector _cVectorRadiusV;
	
	/**
	 * Another vector, orthogonal to the vector along the axis, used for
	 * converting back to Cartesian coordinates.
	 */
	private ContinuousVector _cVectorRadiusW;
	
	/**
	 * The radius of this cylinder. Equivalent to _cVectorRadius.norm().
	 */
	private Double _radius;
	
	/**
	 * The length of this cylinder. Equivalent to _cVectorAlongAxis.norm().
	 */
	private Double _length;
	
	/**
	 * Whether the inside of the cylinder is the inside (true) or the outside
	 * (false) of the domain. 
	 */
	private Boolean _interiorMatchesDomain;
	
	/**
	 * 
	 * @param dPointCenterBase
	 * @param dVectorAlongAxis
	 * @param radius
	 * @param interiorMatchesDomain
	 * @param resolution Grid resolution of the domain.
	 */
	public Cylindrical(DiscreteVector dPointCenterBase, 
			  				DiscreteVector dVectorAlongAxis, Double radius,
			  				Boolean interiorMatchesDomain, Double resolution)
	{
		_dPointCenterBase = new DiscreteVector(dPointCenterBase);
		_dVectorAlongAxis = new DiscreteVector(dVectorAlongAxis);
		_radius = radius;
		_interiorMatchesDomain = interiorMatchesDomain;
		init(resolution);
	}
	
	/**
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		_dPointCenterBase = new DiscreteVector(shapeRoot.getParamParser("pointCenter"));
		_dVectorAlongAxis = new DiscreteVector(shapeRoot.getParamParser("vectorAxis"));
		_radius = shapeRoot.getParamLength("radius");
		_interiorMatchesDomain = shapeRoot.getParamBool("interiorMatchesDomain");
		
		Double res = aDomain.getResolution();
		init(res);
	}
	
	private void init(Double res)
	{
		_dVectorRadiusV = new DiscreteVector();
		_dVectorRadiusW = new DiscreteVector();
		_dVectorAlongAxis.orthoVector(_dVectorRadiusV, _dVectorRadiusW);
		
		_cPointCenterBase = new ContinuousVector();
		_cPointCenterBase.set(_dPointCenterBase, res);
		
		_cVectorAlongAxis = new ContinuousVector();
		_cVectorAlongAxis.set(_dVectorAlongAxis, res);
		
		_length = _cVectorAlongAxis.norm();
		
		_cVectorRadiusV = new ContinuousVector();
		_cVectorRadiusV.set(_dVectorRadiusV);
		_cVectorRadiusV.normalizeVector(_radius);
		
		_cVectorRadiusW = new ContinuousVector();
		_cVectorRadiusW.set(_dVectorRadiusW);
		_cVectorRadiusW.normalizeVector(_radius);
		
		/* IMPORTANT:
		 * - if voronoiPrimary = azimuthCoord, _maxStar = 2.0*Math.PI
		 * - if voronoiPrimary = zenithCoord, _maxStar = Math.PI
		 */
		_voronoiPrimary = azimuthCoord;
		_voronoiSecondary = heightCoord;
		_minPrimary = 0.0;
		_maxPrimary = 2.0*Math.PI;
	}
	
	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		Double[] position = convertToLocal(point);
		Boolean isInsideCylinder = (position[heightCoord] >= 0.0) &&
									(position[heightCoord] <= _length) && 
									(position[radialCoord] <= _radius);
		return Boolean.logicalXor(isInsideCylinder, _interiorMatchesDomain);
	}
	
	public DiscreteVector getRelativePosition(DiscreteVector coord)
	{
		DiscreteVector pointOnPlaneToPoint = new DiscreteVector();
		pointOnPlaneToPoint.sendDiff(coord, _dPointCenterBase);
		return pointOnPlaneToPoint;
	}
	
	public DiscreteVector getAbsolutePosition(DiscreteVector coord)
	{
		DiscreteVector out = new DiscreteVector();
		out.sendSum(coord, _dPointCenterBase);
		return out;
	}
	
	public ContinuousVector getRelativePosition(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.subtract(_cPointCenterBase);
		return out;
	}
	
	public ContinuousVector getAbsolutePosition(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.add(_cPointCenterBase);
		return out;
	}
	
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		Double[] p = convertToLocal(ccIn);
		p[heightCoord] = Math.max(Math.min(p[heightCoord], _length), 0.0);
		p[radialCoord] = _radius;
		// The angle, p[azimuthCoord], doesn't change.
		convertToCartesian(p, ccOut);
	}
	
	public void orthoProj(DiscreteVector dcIn, DiscreteVector dcOut)
	{
		// TODO
	}
	
	/**
	 * 
	 * @param point 
	 * @return 3 Doubles: (0) height, (1) radius, (2) angle
	 */
	public Double[] convertToLocal(ContinuousVector point)
	{
		Double[] out = new Double[3];
		ContinuousVector baseToPoint = getRelativePosition(point);
		Double cosAngle = _cVectorAlongAxis.cosAngle(baseToPoint);
		Double dist = baseToPoint.norm();
		if ( cosAngle.equals(0.0) )
			return null; //TODO
		// cosAngle = dist/height
		out[heightCoord] = dist/cosAngle;
		out[radialCoord] = ExtraMath.triangleSide(dist, out[0]);
		// Now find the position on the axis that is closest to the point.
		ContinuousVector nearestOnAxis = new ContinuousVector(_cVectorAlongAxis);
		nearestOnAxis.times(out[heightCoord]/_length);
		/* 
		 * By subtracting this from the baseToPoint vector, we flatten it onto
		 *  the base circle and find the angle with the reference vector.
		 */
		baseToPoint.subtract(nearestOnAxis);
		out[azimuthCoord] = _cVectorRadiusV.angle(baseToPoint);
		return out;
	}
	
	/**
	 * 
	 * @param coords
	 * @param out
	 */
	private void convertToCartesian(Double[] coords, ContinuousVector out)
	{
		out.set(_cPointCenterBase);
		ContinuousVector temp = new ContinuousVector();
		
		temp.set(_cVectorAlongAxis);
		temp.times(coords[heightCoord]/_length);
		out.add(temp);
		
		temp.set(_cVectorRadiusV);
		temp.times(coords[radialCoord]*Math.cos(coords[azimuthCoord]));
		out.add(temp);
		
		temp.set(_cVectorRadiusW);
		temp.times(coords[radialCoord]*Math.sin(coords[azimuthCoord]));
		out.add(temp);
	}
	
	/**
	 * 
	 * @param coords
	 * @return
	 */
	private ContinuousVector convertToCartesian(Double[] coords)
	{
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(coords, out);
		return out;
	}
	
	/**
	 * 
	 * 
	 * Assumes both points on the surface!
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	@Override
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		Double[] p1 = convertToLocal(point1);
		Double[] p2 = convertToLocal(point2);
		
		Double angle = Math.abs(p1[azimuthCoord] - p2[azimuthCoord]);
		angle = Math.min(angle, 2*Math.PI - angle);
		
		return Math.hypot(p1[heightCoord] - p2[heightCoord], _radius*(angle));
	}
	
	
	@Override
	public Edge bisect(Site site1, Site site2)
	{
		
		return null;
	}
	
	public Vertex intersect(HalfEdge he1, HalfEdge he2)
	{
		Vertex out = new Vertex();
		// TODO
		return out;
	}
	
	public LinkedList<ContinuousVector> getIntersections(
						ContinuousVector position, ContinuousVector vector)
	{
		LinkedList<ContinuousVector> out = new LinkedList<ContinuousVector>();
		if ( _cVectorAlongAxis.isParallel(vector) )
			return out;
		/* 
		 * Flatten both the relative position and the vector to the plane of
		 * the base of the cylinder.
		 */
		ContinuousVector temp = getRelativePosition(position);
		Double[] p = convertToLocal(temp);
		p[heightCoord] = 0.0;
		ContinuousVector pFlat = convertToVector(p);
		Double[] v = convertToLocal(vector);
		v[heightCoord] = 0.0;
		ContinuousVector vFlat = convertToVector(v);
		// Find roots (k) of the equation: |pFlat + k*vFlat| = _radius
		Complex[] roots = ExtraMath.rootsQuadratic(vFlat.prodScalar(vFlat),
							2*vFlat.prodScalar(pFlat),
							pFlat.prodScalar(pFlat) - ExtraMath.sq(_radius));
		// Check roots are real, i.e. the line does intersect the cylinder.
		if ( roots[0].isReal() )
		{
			for ( Complex root : roots )
			{
				// Find the (relative) position of the intersection.
				temp.set(vector);
				temp.times(root.getReal());
				temp.add(position);
				/*
				 * Check the intersection is along the (finite) length of the
				 * cylinder.
				 */
				p = convertToLocal(temp);
				if ( p[heightCoord] >= 0 && p[heightCoord] <= _length )
					out.add(temp);
			}
			/*
			 * Check intersections are different (i.e. if the line is
			 * tangential to the cylinder there is only one intersection).
			 */
			if ( out.get(0).equals(out.get(1)) )
				out.remove(1);
		}
		return out;
	}
	
	public void readyToFollowBoundary(SpatialGrid aSG)
	{
		
	}
	
	/**
	 * 
	 */
	public Boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut,
															SpatialGrid aSG)
	{
		return null;
	}
	
	/**
	 * 
	 */
	public ContinuousVector getNormalInside(ContinuousVector cc)
	{
		return null;
	}
	
	/**
	 * 
	 */
	public Double getDistance(IsShape aBoundary)
	{
		return null;
	}
	
	/**
	 * 
	 */
	public DiscreteVector getNormalDiscrete()
	{
		return null;
	}

	@Override
	public ContinuousVector convertToVector(Double[] local)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Boolean isOutside(DiscreteVector coord)
	{
		// TODO Auto-generated method stub
		return null;
	}
}
