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
			if ( coefficient[1] != 1.0 )
				out += coefficient[1].toString()+"*";
			out += "v ";
		}
		out += "= "+coefficient[2].toString();
		return out;
	}
}
