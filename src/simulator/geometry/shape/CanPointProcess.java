package simulator.geometry.shape;

import simulator.geometry.ContinuousVector;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Vertex;

public interface CanPointProcess extends IsShape
{
	public Double distance(ContinuousVector point1, ContinuousVector point2);
	
	public Edge bisect(Site site1, Site site2);
	
	public Vertex intersect(HalfEdge he1, HalfEdge he2);
	
	public int compare(ContinuousVector point1,
										ContinuousVector point2);
	
}