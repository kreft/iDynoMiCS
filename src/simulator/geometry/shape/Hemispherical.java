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

public class Hemispherical extends IsShape
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * In local coordinates, this gives the index of the azimuthal angle
	 * (if the apex of the hemisphere is equivalent to the North Pole, this is
	 * equivalent of the longitude).
	 */
	public final static int azimuthCoord = 0;
	
	/**
	 * In local coordinates, this gives the index of the azimuth angle, also
	 * known as the polar angle or colatitude (if the apex of the hemisphere
	 * is equivalent to the North Pole, this is equivalent of 90 degrees minus
	 * the latitude).
	 */
	public final static int zenithCoord = 1;
	
	/**
	 * In local coordinates, this gives the index of the radius (distance from
	 * the centre).
	 */
	public final static int radialCoord = 2;
	
	/**
	 * A point at one end of the the cylinder axis.
	 */
	private DiscreteVector _dPointCenterBase;
	
	/**
	 * A vector from _dPointCenterBase to the other end of the cylinder axis.
	 */
	private DiscreteVector _dVectorToApex;
	
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
	
	/**
	 * A point at one end of the the cylinder axis.
	 */
	private ContinuousVector _cPointCenterBase;
	
	/**
	 * A vector from _dPointCenterBase to the other end of the cylinder axis.
	 */
	private ContinuousVector _cVectorToApex;
	
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
	 * Whether the inside of the hemisphere is the inside (true) or the outside
	 * (false) of the domain. 
	 */
	private Boolean _interiorMatchesDomain;
	
	/**
	 * Generic constructor for instance creator (leave empty).
	 */
	public Hemispherical()
	{
		
	}
	
	/**
	 * 
	 * @param dPointCenterBase
	 * @param dVectorToApex
	 * @param interiorMatchesDomain
	 */
	public Hemispherical(DiscreteVector dPointCenterBase,
				DiscreteVector dVectorToApex, Boolean interiorMatchesDomain)
	{
		_dPointCenterBase = new DiscreteVector(dPointCenterBase);
		_dVectorToApex = new DiscreteVector(dVectorToApex);
		_interiorMatchesDomain = interiorMatchesDomain;
	}
	
	/**
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		_dPointCenterBase =
				new DiscreteVector(shapeRoot.getParamParser("pointCenter"));
		_dVectorToApex = 
				new DiscreteVector(shapeRoot.getParamParser("vectorAxis"));
		_interiorMatchesDomain = 
							shapeRoot.getParamBool("interiorMatchesDomain");
		init(aDomain.getResolution());
	}
	
	private void init(Double res)
	{
		_dVectorRadiusV = new DiscreteVector();
		_dVectorRadiusW = new DiscreteVector();
		_dVectorToApex.orthoVector(_dVectorRadiusV, _dVectorRadiusW);
		
		_cPointCenterBase = new ContinuousVector();
		_cPointCenterBase.set(_dPointCenterBase, res);
		
		_cVectorToApex = new ContinuousVector();
		_cVectorToApex.set(_dVectorToApex, res);
		
		_radius = _cVectorToApex.norm();
		
		_cVectorRadiusV = new ContinuousVector();
		_cVectorRadiusV.set(_dVectorRadiusV);
		_cVectorRadiusV.normalizeVector(_radius);
		
		_cVectorRadiusW = new ContinuousVector();
		_cVectorRadiusV.set(_dVectorRadiusW);
		_cVectorRadiusV.normalizeVector(_radius);
		
		/* IMPORTANT:
		 * - if voronoiPrimary = azimuthCoord, _maxStar = Math.PI
		 * - if voronoiPrimary = zenithCoord, _maxStar = 2.0*Math.PI
		 */
		_voronoiPrimary = azimuthCoord;
		_voronoiSecondary = zenithCoord;
		_minPrimary = 0.0;
		_maxPrimary = Math.PI;
	}
	
	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		ContinuousVector baseToPoint = getRelativePosition(point);
		// Work out if the point is inside the hemisphere.
		Boolean isInsideHS = ( _cVectorToApex.cosAngle(baseToPoint) > 0 )
										&& ( baseToPoint.norm() <= _radius);
		// Use the exclusive or (Xor).
		return Boolean.logicalXor(isInsideHS, _interiorMatchesDomain);
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
	
	/**
	 * 
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		ContinuousVector diff = getRelativePosition(ccIn);
		if ( diff.isZero() )
			return;
		Double minDist = Double.MAX_VALUE;
		Double dist;
		for ( ContinuousVector intersection : getIntersections(ccIn, diff) )
		{
			dist = intersection.distance(ccIn);
			if ( dist < minDist )
			{
				ccOut = intersection;
				minDist = dist;
			}
		}
	}
	
	/**
	 * \brief Gets the (Euclidean) distance from given point to the nearest
	 * position on this shape. 
	 */
	public Double getDistance(ContinuousVector point)
	{
		ContinuousVector baseToPoint = getRelativePosition(point);
		Double cosAngle = _cVectorToApex.cosAngle(baseToPoint);
		if ( cosAngle >= 0 )
			return Math.abs(baseToPoint.norm() - _radius);
		/* If the point is below  the hemisphere, use the cosine rule to find
		 * the distance from the circle at the base. Note that, since we want
		 * cos(cosAngle - pi/2) we need to convert it first. 
		 */
		Double angle = Math.abs(Math.acos(cosAngle));
		angle -= Math.PI/2;
		cosAngle = Math.cos(angle);
		Double distSq = baseToPoint.prodScalar(baseToPoint);
		Double out = distSq + ExtraMath.sq(_radius);
		out  += 2 * Math.sqrt(distSq) * cosAngle;
		return Math.sqrt(out);
	}
	
	/**
	 * 
	 * 
	 * 
	 * @param point ContinuousVector describing a point in space (Cartesian
	 * coordinates, relative to global origin) 
	 * @see azimuthCoord, radialCoord, zenithCoord
	 */
	public Double[] convertToLocal(ContinuousVector point)
	{
		return convertToPolar(point, _cVectorRadiusV, _cVectorToApex);
	}
	
	/**
	 * 
	 * 
	 * 
	 * @param point
	 * @param radialVector
	 * @param apexVector
	 * @return
	 * @see azimuthCoord, radialCoord, zenithCoord
	 */
	private Double[] convertToPolar(ContinuousVector point,
					ContinuousVector radialVector, ContinuousVector apexVector)
	{
		Double[] local = new Double[3];
		ContinuousVector baseToPoint = getRelativePosition(point);
		Double height = apexVector.cosAngle(baseToPoint);
		local[radialCoord] = baseToPoint.norm();
		
		// Should be equivalent to Planar.orthoProj() for the base of the HS
		ContinuousVector pointOnPlane = new ContinuousVector(apexVector);
		pointOnPlane.times(-height*local[radialCoord]);
		pointOnPlane.add(baseToPoint);
		
		local[azimuthCoord] = radialVector.angle(pointOnPlane);
		local[zenithCoord] = Math.acos(height);
		return local;
	}
	
	private void convertToCartesian(Double[] local, ContinuousVector out,
							ContinuousVector center, ContinuousVector apexU,
							ContinuousVector radialV, ContinuousVector radialW)
	{
		out.set(center);
		ContinuousVector temp = new ContinuousVector();
		Double k = local[radialCoord] * Math.sin(local[zenithCoord]);
		
		temp.set(radialV);
		temp.times(k * Math.cos(local[azimuthCoord]));
		out.add(temp);
		
		temp.set(radialW);
		temp.times(k * Math.sin(local[azimuthCoord]));
		out.add(temp);
		
		temp.set(apexU);
		temp.times(local[radialCoord] * Math.cos(local[zenithCoord]));
		out.add(temp);
	}
	
	/**
	 * 
	 * 
	 * @param local
	 * @param center
	 * @param apexU
	 * @param radialV
	 * @param radialW
	 * @return
	 */
	private ContinuousVector convertToCartesian(Double[] local,
							ContinuousVector center, ContinuousVector apexU,
							ContinuousVector radialV, ContinuousVector radialW)
	{
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(local, out, center, apexU, radialV, radialW);
		return out;
	}
	
	/**
	 * \brief 
	 * 
	 * @param local Double array describing a position in local coordinates
	 * (spherical, relative to _cPointCenterBase)
	 * @return ContinuousVector describing the position in global coordinates
	 * (Cartesian, relative to global origin)
	 */
	@Override
	public ContinuousVector convertToVector(Double[] local)
	{
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(local, out, _cPointCenterBase, _cVectorToApex,
											_cVectorRadiusV, _cVectorRadiusW);
		return out;
	}
	
	/**
	 * 
	 * Assumes both points on the surface!
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		return _radius * getRelativePosition(point1).angle(getRelativePosition(point2));
	}
	
	/**
	 * 
	 * @param point Any point in Cartesian space.
	 * @return	The vector from _cPointCenterBase to this point.
	 */
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
	
	private ContinuousVector midpoint(Site site1, Site site2)
	{
		ContinuousVector baseToS1 = getRelativePosition(site1);
		ContinuousVector baseToS2 = getRelativePosition(site2);
		ContinuousVector orthog = baseToS1.crossProduct(baseToS2);
		
		Double[] s2 = convertToPolar(site2, baseToS1, orthog);
		s2[0] *= 0.5;
		
		ContinuousVector newW = baseToS1.crossProduct(orthog);
		newW.turnAround();
		newW.normalizeVector(_radius);
		
		return convertToCartesian(s2, _cPointCenterBase, orthog, baseToS1, newW);
	}
	
	/**
	 * 
	 */
	public Edge bisect(Site site1, Site site2)
	{
		// If the sites are co-localised, there's no unique bisector.
		if ( site1.equals(site2) )
			return null;
		Edge out = new Edge();
		// Set the regions on either side of the Edge.
		out.site[0] = site1;
		out.site[1] = site2;
		// Already implicitly the case, but stated explicitly for clarity.
		out.endPoint[0] = null;
		out.endPoint[1] = null;
		
		ContinuousVector baseToS1 = getRelativePosition(site1);
		ContinuousVector baseToS2 = getRelativePosition(site2);
		ContinuousVector orthog = baseToS1.crossProduct(baseToS2);
		orthog.normalizeVector(_radius);
		
		Double[] s = convertToPolar(site2, baseToS1, orthog);
		s[0] *= 0.5;
		
		ContinuousVector midpoint = new ContinuousVector(); 
		convertToCartesian(s, midpoint, _cPointCenterBase, orthog, baseToS1, baseToS1.crossProduct(orthog));
		// TODO
		return out;
	}
	
	public Vertex intersect(Edge edge1, Edge edge2)
	{
		Vertex out = new Vertex();
		// TODO
		return out;
	}
	
	public LinkedList<ContinuousVector> getIntersections(ContinuousVector position,
													ContinuousVector vector) 
	{
		ContinuousVector baseToP = getRelativePosition(position);
		ContinuousVector temp = new ContinuousVector();
		LinkedList<ContinuousVector> out = new LinkedList<ContinuousVector>();
		// Find roots (k) of the equation: |baseToP + k*vector| = _radius
		Complex[] roots = ExtraMath.rootsQuadratic(vector.prodScalar(vector),
						2*vector.prodScalar(baseToP),
						baseToP.prodScalar(baseToP) - ExtraMath.sq(_radius));
		// Check roots are real, i.e. the line does intersect the (hemi)sphere.
		if ( roots[0].isReal() )
		{
			for ( Complex root : roots )
			{
				// Find the (relative) position of the intersection.
				temp.set(vector);
				temp.times(root.getReal());
				temp.add(baseToP);
				// Check this is on the hemisphere, not just the sphere.
				if ( _cVectorToApex.cosAngle(temp) >= 0.0 )
				{
					temp.add(_cPointCenterBase);
					out.add(new ContinuousVector(temp));
				}
			}
			/*
			 * Check intersections are different (i.e. line not tangential to
			 * the hemisphere).
			 */
			if ( out.get(0).equals(out.get(1)) )
				out.remove(1);
		}
		return out;
	}
	
	public void readyToFollowBoundary(SpatialGrid aSG)
	{
		// TODO Auto-generated method stub
	}
	
	public Boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut,
															SpatialGrid aSG) 
	{
		return false;
	}
	
	
	public ContinuousVector getNormalInside()
	{
		return null;
	}
	
	public Double getDistance(IsShape aBoundary)
	{
		
		return null;
	}
	
	public DiscreteVector getNormalDiscrete()
	{
		return null;
	}

	@Override
	public void orthoProj(DiscreteVector dcIn, DiscreteVector dcOut)
	{
		// TODO Auto-generated method stub
	}

	

	@Override
	public Boolean isOutside(DiscreteVector coord)
	{
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public ContinuousVector getEdgePointFromPrimary(Edge edge, Double primaryValue)
	{
		// TODO Auto-generated method stub
		return null;
	}
	
	public StringBuffer writeShapeInformation(StringBuffer outputString)
	{
		outputString.append("<Surface shape=\"Hemispherical\"");
		outputString.append(" cPointCenterBase=\""+
										_cPointCenterBase.toString()+"\"");
		outputString.append(" cVectorRadiusV=\""+
											_cVectorRadiusV.toString()+"\"");
		outputString.append(" cVectorRadiusW=\""+
											_cVectorRadiusW.toString()+"\"");
		outputString.append("/>\n");
		return outputString;
	}
	
	public StringBuffer getSitesHeader()
	{
		return new StringBuffer("azimuth,zenith");
	}
	
	public StringBuffer getEdgesHeader()
	{
		return new StringBuffer("azimuth1,zenith1,azimuth2,zenith2");
	}

	@Override
	public void clipEdgeToLimits(Edge edge) {
		// TODO Auto-generated method stub
		
	}
}
