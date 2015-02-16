package simulator.geometry.pointProcess;

public class HalfEdge
{
	public Edge edge;
	
	public HalfEdge leftNeighbor, rightNeighbor;
	
	public Boolean deleted;
	
	public Vertex vertex;
	
	private static final int left = 0;
	private static final int right = 1;
	
	/**
	 * Integer denoting whether this HalfEdge is on the left (0) or the right
	 * (1) hand side of the Edge it corresponds to. 
	 */
	private int leftRight;
	
	/**
	 * leftRight is set to left by default.
	 */
	public HalfEdge()
	{
		this.edge = null;
		this.vertex = null;
		this.leftRight = left;
		this.deleted = false;
	}
	
	/**
	 * \brief 
	 * 
	 * 
	 * @param edge Edge to be associated with this HalfEdge.
	 * @param placeOnLeft Boolean noting whether to place this HalfEdge on the
	 * left (true) or on the right (false).
	 */
	public HalfEdge(Edge edge, Boolean placeOnLeft)
	{
		this.edge = edge;
		this.vertex = null;
		this.leftRight = placeOnLeft ? left : right;
		this.deleted = false;
	}
	
	public Boolean isOnLeft()
	{
		return leftRight == left;
	}
	
	public Boolean isOnRight()
	{
		return leftRight == right;
	}
	
	public Site getLeftRegion()
	{
		return edge.region[left];
	}
	
	public Site getRightRegion()
	{
		return edge.region[right];
	}
	
	public Boolean isNearVertical()
	{
		return edge.coefficient[0] == 1.0;
	}
	
	public String toString()
	{
		String out = "HalfEdge ";
		if ( deleted )
			out += "(deleted) ";
		
		out += "on the ";
		out += ( isOnLeft() ) ? "left" : "right";
		
		out += " with ";
		if ( leftNeighbor == null )
		{
			if ( rightNeighbor == null )
				out += "no neighbors";
			else
				out += "right neighbor only";
		}
		else if ( rightNeighbor == null )
			out += "left neighbor only";
		else
			out += "left and right neighbors";
		
		out += " and ";
		
		if ( edge == null )
			return out + "no edge attached";
		
		return out + edge.toString();
	}
}
