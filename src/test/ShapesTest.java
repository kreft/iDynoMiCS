package test;

import simulator.geometry.*;
import simulator.geometry.shape.*;

public class ShapesTest
{
	
	public static void main(String[] args) 
	{
		Double resolution = 1.0;
		
		DiscreteVector pointOnPlane = new DiscreteVector();
		DiscreteVector vectorOut = new DiscreteVector(1, 0, 0);
		
		Planar plane1 = new Planar(pointOnPlane, vectorOut, resolution);
		
		
		
	}
	
	
}
