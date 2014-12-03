package simulator.geometry.shape;

import simulator.geometry.ContinuousVector;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.Site;

public interface CanPointProcess extends IsShape
{
	public Double distance(ContinuousVector point1, ContinuousVector point2);
	
	public Edge bisect(Site site1, Site site2);
}