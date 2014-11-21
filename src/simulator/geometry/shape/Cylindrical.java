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
	private DiscreteVector _dPointOnAxis;
	
	
	/**
	 * The radius of this cylinder.
	 */
	private Double _radius;
	
	private Double length;
	
	/**
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		
		
	}
	
	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector cV)
	{
		
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
