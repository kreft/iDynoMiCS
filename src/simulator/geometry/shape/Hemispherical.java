package simulator.geometry.shape;

import java.io.Serializable;

import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.Domain;
import utils.ExtraMath;
import utils.XMLParser;

public class Hemispherical implements IsShape, Serializable
{
	/**
	 * A point on the cylinder axis
	 */
	private DiscreteVector _dPointCenterBase;
	
	/**
	 * 
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
	
	
	private ContinuousVector _cPointCenterBase;
	
	
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
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		_dPointCenterBase = new DiscreteVector(shapeRoot.getParamParser("pointCenter"));
		_dVectorToApex = new DiscreteVector(shapeRoot.getParamParser("vectorAxis"));
		
		_dVectorRadiusV = new DiscreteVector();
		_dVectorRadiusW = new DiscreteVector();
		_dVectorToApex.orthoVector(_dVectorRadiusV, _dVectorRadiusW);
		
		Double res = aDomain.getResolution();
		_cPointCenterBase = new ContinuousVector(_dPointCenterBase, res);
		_cVectorToApex = new ContinuousVector(_dVectorToApex, res);
		_radius = _cVectorToApex.norm();
		_cVectorRadiusV = new ContinuousVector(_dVectorRadiusV, res);
		_cVectorRadiusV.normalizeVector(_radius);
		_cVectorRadiusW = new ContinuousVector(_dVectorRadiusW, res);
		_cVectorRadiusV.normalizeVector(_radius);
		
		_interiorMatchesDomain = shapeRoot.getParamBool("interiorMatchesDomain");
	}

	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		ContinuousVector baseToPoint = baseToPoint(point);
		// Work out if the point is inside the hemisphere.
		Boolean insideHS = ( _cVectorToApex.cosAngle(baseToPoint) > 0 ) 
									&& ( baseToPoint.norm() <= _radius);
		return ( _interiorMatchesDomain )? insideHS : !insideHS;
	}
	
	/**
	 * 
	 */
	public ContinuousVector intersection(ContinuousVector position,
													ContinuousVector vector)
	{
		
		return null;
	}
	
	/**
	 * 
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		
		
	}
	
	/**
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		
		return null;
	}
	
	/**
	 * 
	 */
	public Double getDistance(ContinuousVector point)
	{
		ContinuousVector baseToPoint = baseToPoint(point);
		Double cosAngle = _cVectorToApex.cosAngle(baseToPoint);
		if ( cosAngle >= 0 )
			return baseToPoint.norm() - _radius;
		/* If the point is below  the hemisphere, use the cosine rule to find
		 * the distance from the circle at the base. Note that, since we want
		 * cos(cosAngle - pi/2) we need to convert it first. 
		 */
		Double angle = Math.abs(Math.acos(cosAngle));
		angle -= Math.PI/2;
		cosAngle = Math.cos(angle);
		Double distSq = baseToPoint.normSq();
		Double out = distSq + ExtraMath.sq(_radius);
		out  += 2 * Math.sqrt(distSq) * cosAngle;
		return Math.sqrt(out);
	}
	
	/**
	 * 
	 * @param point Must be on (round) surface of hemisphere!
	 * @return 2 Doubles: (0) angle on the circular base, (1) angle with the
	 * vector to the apex.
	 */
	private Double[] convertToPolar(ContinuousVector point)
	{
		Double[] out = new Double[2];
		
		
		out[1] = _cVectorToApex.cosAngle(point);
		
		ContinuousVector pointOnPlane = new ContinuousVector(_cVectorToApex);
		pointOnPlane.times(-out[1]);
		pointOnPlane.add(point);
		
		out[0] = _cVectorRadiusV.angle(pointOnPlane);
		
		out[1] = Math.acos(out[1]);
		
		return out;
	}
	
	/**
	 * 
	 * @param thetaPhi
	 * @param out
	 */
	private void convertToCartesian(Double[] thetaPhi, ContinuousVector out)
	{
		out.set(_cPointCenterBase);
		ContinuousVector temp = new ContinuousVector(_cVectorRadiusV);
		temp.times(Math.cos(thetaPhi[0]));
		out.add(temp);
		temp.set(_cVectorToApex);
		temp.times(Math.cos(thetaPhi[1]));
		out.add(temp);
	}
	
	/**
	 * 
	 * @param heightTheta
	 * @return
	 */
	private ContinuousVector convertToCartesian(Double[] heightTheta)
	{
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(heightTheta, out);
		return out;
	}
	
	/**
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		Double[] p1 = convertToPolar(point1);
		Double[] p2 = convertToPolar(point2);
		// First check that the angle around the circular base is correct.
		Double out = Math.abs(p1[0] - p2[0]);
		out = Math.min(out, 2*Math.PI - out);
		/* Then calculate the Great-Circle Distance using the  
		 * Spherical Law of Cosines.
		 */
		out = Math.sin(out);
		out *= Math.cos(p1[1]) * Math.cos(p2[1]);
		out += Math.sin(p1[1]) * Math.sin(p2[1]);
		return _radius * Math.acos(out);
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
}
