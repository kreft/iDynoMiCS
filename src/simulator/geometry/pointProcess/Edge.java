package simulator.geometry.pointProcess;

public class Edge
{
	
	public Vertex[] endPoint;
	
	
	public Site[] region;
	
	public Double[] coefficient;
	
	public Edge()
	{
		endPoint = new Vertex[2];
		region = new Site[2];
		coefficient = new Double[3];
	}
	
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
		out += "= "+coefficient[2].toString()+" [L: ";
		if ( region[0] == null )
			out += "empty";
		else
			out += region[0].toString();
		out += ", R: ";
		if ( region[1] == null )
			out += "empty";
		else
			out += region[1].toString();
		
		return out + "]";
	}
}
