package test;

import simulator.geometry.shape.*;
import utils.XMLParser;

public class VoronoiTest
{
	public static void main(String[] args) 
	{
		XMLParser protocolFile = new XMLParser(args[0]);
		
		XMLParser shapeRoot = protocolFile.getChildParser("shape");
		
		
		
	}
}
