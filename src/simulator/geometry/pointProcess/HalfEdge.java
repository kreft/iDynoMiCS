package simulator.geometry.pointProcess;

public class HalfEdge
{
	public Edge edge;
	
	public HalfEdge leftNeighbor, rightNeighbor;
	
	public Boolean deleted;
	
	public Vertex vertex;
	
	public static final int left = 0;
	public static final int right = 1;
	
	/**
	 * Integer denoting whether this HalfEdge is on the left (0) or the right
	 * (1) hand side of the Edge it corresponds to. 
	 */
	private int side;
	
	
	/**
	 * side is set to left by default.
	 */
	public HalfEdge()
	{
		this.edge = null;
		this.vertex = null;
		this.side = Edge.left;
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
		this.side = placeOnLeft ? Edge.left : Edge.right;
		this.deleted = false;
	}
	
	public int getSide()
	{
		return side;
	}
	
	public Boolean isOnLeft()
	{
		return side == Edge.left;
	}
	
	public Boolean isOnRight()
	{
		return side == Edge.right;
	}
	
	/**
	 * \brief Gets the site below this HalfEdge's associated Edge.
	 * 
	 * @return Site below the Edge associated with this HalfEdge. Null if this
	 * Edge is not yet set.
	 */
	public Site getSiteBelow()
	{
		return ( this.edge == null ) ? null : this.edge.getSiteBelow();
	}
	
	/**
	 * \brief Gets the site above this HalfEdge's associated Edge.
	 * 
	 * @return Site above the Edge associated with this HalfEdge. Null if this
	 * Edge is not yet set.
	 */
	public Site getSiteAbove()
	{
		return ( this.edge == null ) ? null : this.edge.getSiteAbove();
	}
	
	/**
	 * 
	 * 
	 * @return
	 */
	public Site getSiteOnLeft()
	{
		return this.isOnLeft() ? this.getSiteBelow() : this.getSiteAbove();
	}
	
	/**
	 * 
	 * 
	 * @return
	 */
	public Site getSiteOnRight()
	{
		return this.isOnRight() ? this.getSiteBelow() : this.getSiteAbove();
	}
	
	public Boolean isNearVertical()
	{
		return edge.isNearVertical();
	}
	
	/**
	 * TODO Maybe move to space?
	 * 
	 * @param primaryValue
	 * @return
	 */
	public Double getSecondaryValue(Double primaryValue)
	{
		return edge.coefficient[2] - edge.coefficient[0]*primaryValue;
	}
	
	/**
	 * \brief Prints out a description of this HalfEdge to screen.
	 * 
	 * Useful during testing and debugging.
	 */
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
