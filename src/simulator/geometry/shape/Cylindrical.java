package simulator.geometry.shape;

import java.io.Serializable;

import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.Domain;
import utils.XMLParser;

public class Cylindrical implements IsShape, Serializable
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
	
	
	private ContinuousVector _cPointCenterBase;
	
	
	private ContinuousVector _cVectorAlongAxis;
	
	/**
	 * The radius of this cylinder.
	 */
	private Double _radius;
	
	/**
	 * The length of this cylinder. Equivalent to _cVectorAlongAxis.norm().
	 */
	private Double _length;
	
	/**
	 * Whether the inside of the cyclinder is the inside (true) or the outside
	 * (false) of the domain. 
	 */
	private Boolean _interiorMatchesDomain;
	
	/**
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		_dPointCenterBase = new DiscreteVector(shapeRoot.getParamParser("pointCenter"));
		_dVectorAlongAxis = new DiscreteVector(shapeRoot.getParamParser("vectorAxis"));
		_radius = shapeRoot.getParamLength("radius");
		
		Double res = aDomain.getResolution();
		_cPointCenterBase = new ContinuousVector(_dPointCenterBase, res);
		_cVectorAlongAxis = new ContinuousVector(_dVectorAlongAxis, res);
		_length = _cVectorAlongAxis.norm();
		
		_interiorMatchesDomain = shapeRoot.getParamBool("interiorMatchesDomain");
	}
	
	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		ContinuousVector temp = new ContinuousVector(point);
		
		
		return null;
	}
	
	/**
	 * 
	 */
	public Boolean isOnBoundary(ContinuousVector cV, Double res)
	{
		
		return null;
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
	public Double getDistance(IsShape aShape)
	{
		
		return null;
	}
	
	/**
	 * 
	 */
	public Double getDistance(ContinuousVector cc)
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
	public void readyToFollowBoundary(SpatialGrid aSG)
	{
		
		
	}
	
	/**
	 * 
	 */
	public Boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut,
															SpatialGrid aSG)
	{
		
		return false;
	}
	
	/**
	 * 
	 */
	public DiscreteVector getNormalDC()
	{
		
		return null;
	}
	
	
}
