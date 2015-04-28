package simulator.geometry.pointProcess;

public class HalfEdge
{
	public Edge edge;
	
	public HalfEdge previousNeighbor, nextNeighbor;
	
	public Boolean deleted;
	
	/**
	 * Boolean denoting whether this HalfEdge is on the out-bound (true) or
	 * the in-bound (false) side of the Edge it corresponds to. 
	 */
	private Boolean isOutbound;
	
	/**
	 * side is set to left by default.
	 */
	public HalfEdge(Boolean isOutbound)
	{
		this.edge = null;
		this.isOutbound = isOutbound;
		this.deleted = false;
	}
	
	/**
	 * \brief 
	 * 
	 * 
	 * @param edge Edge to be associated with this HalfEdge.
	 * @param isOutbound Boolean noting whether to make this HalfEdge 
	 * out-bound (true) or in-bound (false).
	 */
	public HalfEdge(Edge edge, Boolean isOutbound)
	{
		this.edge = edge;
		this.isOutbound = isOutbound;
		this.deleted = false;
	}
	
	public Boolean isOutbound()
	{
		return this.isOutbound;
	}
	
	public Boolean isInbound()
	{
		return ! this.isOutbound;
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
	public Site getSiteBehind()
	{
		return this.isOutbound ? this.getSiteBelow() : this.getSiteAbove();
	}
	
	/**
	 * 
	 * 
	 * @return
	 */
	public Site getSiteAhead()
	{
		return this.isOutbound ? this.getSiteAbove() : this.getSiteBelow();
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
		
		out += ( isOutbound ) ? "outbound" : "inbound";
		
		out += " with ";
		if ( previousNeighbor == null )
		{
			if ( nextNeighbor == null )
				out += "no neighbors";
			else
				out += "next neighbor only";
		}
		else if ( nextNeighbor == null )
			out += "previous neighbor only";
		else
			out += "both neighbors";
		
		out += " and ";
		
		if ( edge == null )
			return out + "no edge attached";
		
		return out + edge.toString();
	}
}
