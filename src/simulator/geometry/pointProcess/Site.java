package simulator.geometry.pointProcess;

import simulator.geometry.ContinuousVector;

public class Site extends ContinuousVector
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public Site(ContinuousVector cV)
	{
		super(cV);
	}
	
	public Site(Double x, Double y, Double z)
	{
		super(x, y, z);
	}
}