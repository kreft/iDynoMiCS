package test;

import java.util.LinkedList;

import simulator.geometry.DiscreteVector;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.shape.Planar;

public class VoronoiTest
{
	public static void main(String[] args) 
	{
		Double resolution = 1.0;
		
		DiscreteVector origin = new DiscreteVector(0, 0, 0);
		DiscreteVector tenTen = new DiscreteVector(0, 25, 25);
		DiscreteVector iUnit = new DiscreteVector(1, 0, 0);
		DiscreteVector jUnit = new DiscreteVector(0, 1, 0);
		DiscreteVector kUnit = new DiscreteVector(0, 0, 1);
		
		Planar mainShape = new Planar(origin, iUnit, resolution);
		System.out.println("This "+mainShape.toString()+"\n");
		
		LinkedList<Planar> walls = new LinkedList<Planar>();
		walls.add(new Planar(origin, jUnit, resolution));
		walls.add(new Planar(origin, kUnit, resolution));
		walls.add(new Planar(tenTen, jUnit, resolution));
		walls.add(new Planar(tenTen, kUnit, resolution));
		mainShape.restrictPlane(walls);
		
		LinkedList<Site> sites = new LinkedList<Site>();
		sites.add(new Site(0.0,206.2160890468038,127.96839725139543))
	}
}
