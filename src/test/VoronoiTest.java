package test;

import java.util.LinkedList;

import simulator.geometry.DiscreteVector;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Voronoi;
import simulator.geometry.shape.Planar;

public class VoronoiTest
{
	public static void main(String[] args) 
	{
		Double resolution = 8.0;
		
		DiscreteVector origin = new DiscreteVector(0, 0, 0);
		DiscreteVector tenTen = new DiscreteVector(0, 33, 33);
		DiscreteVector iUnit = new DiscreteVector(1, 0, 0);
		DiscreteVector jUnit = new DiscreteVector(0, 1, 0);
		DiscreteVector kUnit = new DiscreteVector(0, 0, 1);
		
		Planar mainShape = new Planar(origin, iUnit, resolution);
		//System.out.println("This "+mainShape.toString()+"\n");
		
		LinkedList<Planar> walls = new LinkedList<Planar>();
		walls.add(new Planar(origin, jUnit, resolution));
		walls.add(new Planar(origin, kUnit, resolution));
		walls.add(new Planar(tenTen, jUnit, resolution));
		walls.add(new Planar(tenTen, kUnit, resolution));
		mainShape.restrictPlane(walls);
		
		LinkedList<Site> sites = new LinkedList<Site>();
		/*
		sites.add(new Site(0.0,5.0,8.0));
		sites.add(new Site(0.0,2.0,2.0));
		sites.add(new Site(0.0,6.0,4.0));
		*/
		
		/*
		sites.add(new Site(0.0,8.0,2.0));
		sites.add(new Site(0.0,4.0,4.0));
		sites.add(new Site(0.0,5.0,8.0));
		*/
		
		// Kieran's sites
		//*
		sites.add(new Site(0.0,206.2160890468038,127.96839725139543));
        sites.add(new Site(0.0,82.32194396558708,56.39845559072465));
        sites.add(new Site(0.0,243.45542809645158,50.96943378402841));
        sites.add(new Site(0.0,1.3251214345680484,171.5992947408901));
        sites.add(new Site(0.0,111.61296908358663,68.67900075137223));
        sites.add(new Site(0.0,62.23473496963386,202.16753976573924));
        sites.add(new Site(0.0,62.67681012406526,131.21477136833795));
        sites.add(new Site(0.0,153.22138366950085,70.28696698109235));
        sites.add(new Site(0.0,218.9406070414186,243.70034090543544));
        sites.add(new Site(0.0,122.74295954348423,65.31491274814576));
        sites.add(new Site(0.0,203.88258221137468,171.84703363844423));
        sites.add(new Site(0.0,138.7268577012289,94.02677757643369));
        sites.add(new Site(0.0,85.28252017116063,252.9746284887051));
        sites.add(new Site(0.0,97.07556563528796,98.5329141436925));
        sites.add(new Site(0.0,237.1318113246396,165.5735247757794));
        sites.add(new Site(0.0,190.94077889935267,17.997758178983638));
        sites.add(new Site(0.0,56.5653477270338,235.34520855466639));
        sites.add(new Site(0.0,91.29195509368002,107.33674222593257));
        sites.add(new Site(0.0,85.34218125016955,151.15281661914483));
        sites.add(new Site(0.0,50.176775152154704,203.39704176232215));
        sites.add(new Site(0.0,180.2941531888178,81.12950704316306));
        sites.add(new Site(0.0,40.17916023980316,63.518423127145105));
        sites.add(new Site(0.0,33.50030733144122,257.4369211402401));
        sites.add(new Site(0.0,244.34243200880843,207.57147045344044));
        sites.add(new Site(0.0,261.856334022953,235.4766980966265));
        sites.add(new Site(0.0,103.00285893152787,92.58139657186584));
        sites.add(new Site(0.0,33.454785653891165,246.59793928017893));
        sites.add(new Site(0.0,83.80663276281275,82.83789198266629));
        sites.add(new Site(0.0,131.87560793529698,57.50981580910492));
        sites.add(new Site(0.0,222.41650655175755,190.290678757574));
        sites.add(new Site(0.0,68.5927181961352,84.59952973822497));
        sites.add(new Site(0.0,151.81728949121702,151.5676829854818));
        sites.add(new Site(0.0,187.14356504246678,250.75143779522674));
        sites.add(new Site(0.0,209.04129817533493,204.29868041769666));
        //*/
        Voronoi voronoi = new Voronoi(mainShape, sites);
	}
}
