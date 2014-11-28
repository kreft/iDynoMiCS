package simulator.geometry.shape;

import simulator.geometry.ContinuousVector;

public interface CanPointProcess extends IsShape
{
	public Double distance(ContinuousVector point1, ContinuousVector point2);
}