package simulator.geometry.pointProcess;

import simulator.geometry.ContinuousVector;

public class Vertex extends ContinuousVector
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * 
	 */
	public Double starValue;
	
	public Vertex()
	{
		super();
		starValue = Double.NaN;
	}
	
	public Vertex(ContinuousVector position)
	{
		super(position);
		starValue = Double.NaN;
	}
}
