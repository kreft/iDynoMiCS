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
	
	/**
	 * The HalfEdge associated with this vertex that allows it to be inserted
	 * into the SweepTable. This will typically be a dummy HalfEdge to begin
	 * with, i.e. without an Edge attached.
	 */
	public HalfEdge previousHE;
	
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
