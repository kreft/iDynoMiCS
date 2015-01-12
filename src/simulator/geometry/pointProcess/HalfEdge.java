package simulator.geometry.pointProcess;

public class HalfEdge
{
	public Edge edge;
	
	public HalfEdge leftNeighbor, rightNeighbor;
	
	public Boolean deleted;
	
	/**
	 * Integer donting whether this HalfEdge is on the left (0) or the right
	 * (1) hand side of the Edge it corresponds to. 
	 */
	public int leftRight;
	
	
}
