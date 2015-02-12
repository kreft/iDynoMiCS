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
}
