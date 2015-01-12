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
import utils.XMLParser;

public class Spherical implements IsShape, CanPointProcess, Serializable
{
	/**
	 * Centre of the sphere
	 */
	private DiscreteVector _dPointCenter;
	
	/**
	 * Centre of the sphere
	 */
	private ContinuousVector _cPointCenter;
	
	/**
	 * The radius of the sphere
	 */
	private Double _radius;
	
	@Override
	public Double distance(ContinuousVector point1, ContinuousVector point2) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Edge bisect(Site site1, Site site2) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public Vertex intersect(HalfEdge he1, HalfEdge he2)
	{
		Vertex out = new Vertex();
		// TODO
		return out;
	}
	
	@Override
	public void readShape(XMLParser shapeRoot, Domain aDomain) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public Boolean isOutside(ContinuousVector cV) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public ContinuousVector getOrthoProj(ContinuousVector ccIn) {
		// TODO Auto-generated method stub
		return null;
	}
	
	private ContinuousVector baseToPoint(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.subtract(_cPointCenter);
		return out;
	}
	
	public Double getDistance(ContinuousVector point) 
	{
		ContinuousVector baseToPoint = baseToPoint(point);
		return Math.abs(baseToPoint.norm() - _radius);
	}
	
	/**
	 * TODO 
	 */
	public final int compare(ContinuousVector point1,
											ContinuousVector point2)
	{
		/*
		Double[] p1 = convertToPolar((ContinuousVector) point1);
		Double[] p2 = convertToPolar((ContinuousVector) point2);
		int out = (int) Math.signum(p1[1] - p2[1]);
		if ( out == 0 )
			out = (int) Math.signum(p1[0] - p2[0]);
		return out;
		*/
		return 0;
	}
}
