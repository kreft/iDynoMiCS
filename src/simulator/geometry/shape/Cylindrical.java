package simulator.geometry.shape;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.Domain;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Vertex;
import utils.ExtraMath;
import utils.XMLParser;

public class Cylindrical implements IsShape, CanPointProcess, Serializable
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	
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
	 * @param res	Grid resolution of the domain.
	 */
	public Cylindrical(DiscreteVector dPointCenterBase, 
			  				DiscreteVector dVectorAlongAxis, Double radius,
			  				Boolean interiorMatchesDomain, Double res)
	{
		_dPointCenterBase = new DiscreteVector(dPointCenterBase);
		_dVectorAlongAxis = new DiscreteVector(dVectorAlongAxis);
		_radius = radius;
		_interiorMatchesDomain = interiorMatchesDomain;
		
		init(res);
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
		_cPointCenterBase = new ContinuousVector(_dPointCenterBase, res);
		_cVectorAlongAxis = new ContinuousVector(_dVectorAlongAxis, res);
		_length = _cVectorAlongAxis.norm();
		_cVectorRadiusV = new ContinuousVector(_dVectorRadiusV, res);
		_cVectorRadiusV.normalizeVector(_radius);
		_cVectorRadiusW = new ContinuousVector(_dVectorRadiusW, res);
		_cVectorRadiusV.normalizeVector(_radius);
	}
	
	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		Double[] position = convertToCylindrical(point);
		Boolean isInsideCylinder = (position[0] >= 0.0) &&
					(position[0] <= _length) && (position[1] <= _radius);
		return Boolean.logicalXor(isInsideCylinder, _interiorMatchesDomain);
	}
	
	/**
	 * 
	 */
	public Boolean isOnBoundary(ContinuousVector point, Double res)
	{
		Double[] position = convertToCylindrical(point);
		if ((position[0] < 0.0) || (position[0] > _length))
			return false;
		if ( _interiorMatchesDomain )
			return ( position[1] >= _radius && position[1] <= _radius + res);
		else 
			return ( position[1] <= _radius && position[1] >= _radius - res);
	}
	
	/**
	 * 
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		Double[] coords = convertToCylindrical(ccIn);
		coords[0] = Math.max(Math.min(coords[0], _length), 0.0);
		coords[1] = _radius;
		convertToCartesian(coords, ccOut);
	}
	
	/**
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		ContinuousVector out = new ContinuousVector();
		orthoProj(ccIn, out);
		return out;
	}
	
	/**
	 * 
	 */
	public Double getDistance(ContinuousVector cc)
	{
		Double[] coords = convertToCylindrical(cc);
		Double d = Math.max(coords[0] - _length, - coords[0]);
		Double r = Math.abs(coords[1] - _radius);
		if ( d < 0.0)
			return Math.hypot(d, r);
		else
			return r;
	}
	
	/**
	 * 
	 * @param point 
	 * @return 3 Doubles: (0) height, (1) radius, (2) angle
	 */
	private Double[] convertToCylindrical(ContinuousVector point)
	{
		Double[] out = new Double[3];
		ContinuousVector baseToPoint = baseToPoint(point);
		Double cosAngle = _cVectorAlongAxis.cosAngle(baseToPoint);
		Double dist = baseToPoint.norm();
		if ( cosAngle.equals(0.0) )
			return null; //TODO
		// cosAngle = dist/height
		out[0] = dist/cosAngle;
		out[1] = ExtraMath.triangleSide(dist, out[0]);
		// Now find the position on the axis that is closest to the point.
		ContinuousVector nearestOnAxis = new ContinuousVector(_cVectorAlongAxis);
		nearestOnAxis.times(out[0]/_length);
		/* By subtracting this from the baseToPoint vector, we flatten it onto
		 *  the base circle and find the angle with the reference vector.
		 */
		baseToPoint.subtract(nearestOnAxis);
		out[2] = _cVectorRadiusV.angle(baseToPoint);
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
		
		ContinuousVector temp = new ContinuousVector(_cVectorAlongAxis);
		temp.times(coords[0]/_length);
		out.add(temp);
		
		temp.set(_cVectorRadiusV);
		temp.times(coords[1]*Math.cos(coords[2]));
		out.add(temp);
		
		temp.set(_cVectorRadiusW);
		temp.times(coords[1]*Math.sin(coords[2]));
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
	 * @param point1
	 * @param point2
	 * @return
	 */
	public Double distanceAlongSurface(ContinuousVector point1,
													ContinuousVector point2)
	{
		Double[] p1 = convertToCylindrical(point1);
		Double[] p2 = convertToCylindrical(point2);
		
		Double angle = Math.abs(p1[1] - p2[1]);
		angle = Math.min(angle, 2*Math.PI);
		
		return Math.hypot(p1[0] - p2[0], _radius*(angle));
	}
	
	/**
	 * 
	 * @param point Any point in Cartesian space.
	 * @return	The vector from _cPointCenterBase to this point.
	 */
	private ContinuousVector baseToPoint(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.subtract(_cPointCenterBase);
		return out;
	}
	
	
	@Override
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		
		return null;
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
	
	/**
	 * TODO Check!
	 */
	public void orderLexicographically(ContinuousVector[] list)
	{
		List<ContinuousVector> events = 
				new ArrayList<ContinuousVector>(list.length);
		
		Collections.sort(events, new Comparator<Object>() {
			@Override
			public int compare(Object point1, Object point2)
			{
				Double[] p1 = convertToCylindrical((ContinuousVector) point1);
				Double[] p2 = convertToCylindrical((ContinuousVector) point2);
				int out = (int) Math.signum(p1[1] - p2[1]);
				if ( out == 0 )
					out = (int) Math.signum(p1[0] - p2[0]);
				return out;
			}
		});

		int i = 0;
		for ( ContinuousVector e : events )
		{
			list[i] = e;
			i++;
		}
	}
	
	/**
	 * TODO Check!
	 */
	public final int compare(ContinuousVector point1,
											ContinuousVector point2)
	{		
		Double[] p1 = convertToCylindrical(point1);
		Double[] p2 = convertToCylindrical(point2);
		int out = (int) Math.signum(p1[1] - p2[1]);
		if ( out == 0 )
			out = (int) Math.signum(p1[0] - p2[0]);
		return out;
	}
}
