package simulator.geometry.pointProcess;

public class Edge
{
	public Vertex[] endPoint;
	
	public Site[] site;
	
	public Double[] coefficient;
	
	public static final int left = 0;
	public static final int right = 1;
	
	public static final int below = 0;
	public static final int above = 1;
	
	public Edge()
	{
		endPoint = new Vertex[2];
		site = new Site[2];
		coefficient = new Double[3];
	}
	
	public Boolean areEndPointsSet()
	{
		return (endPoint[left] != null) && (endPoint[right] != null);
	}
	
	/**
	 * \brief Set the 
	 * 
	 * @param vertex
	 * @param setOnLeft
	 */
	public void setEndPoint(Vertex vertex, Boolean setOnLeft)
	{
		endPoint[setOnLeft ? left : right] = vertex;
	}
	
	public Site getSiteBelow()
	{
		return site[below];
	}
	
	public Site getSiteAbove()
	{
		return site[above];
	}
	
	public Boolean isNearVertical()
	{
		return coefficient[0] == 1.0;
	}
	
	/**
	 * 
	 */
	public String toString()
	{
		String out = "Edge: ";
		if ( coefficient[0] != 0.0 )
		{
			if ( coefficient[0] != 1.0 )
				out += coefficient[0].toString()+"*";
			out += "u ";
		}
		if ( coefficient[1] != 0.0 )
		{
			if ( coefficient[1] > 0.0 )
				out += "+ ";
			if ( coefficient[1] != 1.0 )
				out += coefficient[1].toString()+"*";
			out += "v ";
		}
		out += "= "+coefficient[2].toString()+" [S0: ";
		if ( site[0] == null )
			out += "empty";
		else
			out += site[0].toString();
		out += ", S1: ";
		if ( site[1] == null )
			out += "empty";
		else
			out += site[1].toString();
		
		return out + "]";
	}
}
