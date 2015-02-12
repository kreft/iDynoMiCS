package simulator.geometry.pointProcess;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.ListIterator;

import simulator.geometry.ContinuousVector;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.shape.IsShape;
import utils.LogFile;

public class Voronoi
{
	private LinkedList<ContinuousVector> priorityQueue;
	
	private ListIterator<ContinuousVector> pqIterator;
	
	private SweepTable sweepTable;
	
	private IsShape _space;
	
	private IsShape[] boundaries;
	
	private LinkedList<Edge> edges;
	
	private ContinuousVector nextEvent;
	
	public Voronoi(IsShape space, LinkedList<Site> sites)
	{
		_space = space;
		
		// Initialise the Priority Queue
		priorityQueue = new LinkedList<ContinuousVector>();
		for (ContinuousVector site : sites)
		{
			priorityQueue.add(site);
			LogFile.writeLog("Added to Priority Queue: site at "+site.toString());
		}
		Collections.sort(priorityQueue, new Comparator<ContinuousVector>() {
			public final int compare(ContinuousVector p1, ContinuousVector p2)
				{ return _space.compare(p1, p2);}});
		LogFile.writeLog("Priority Queue sorted");
		pqIterator = priorityQueue.listIterator();
		
		// Initialise the Sweep Table
		sweepTable = new SweepTable(_space, sites.size());
		
		HalfEdge leftBoundary, rightBoundary, bisector;
		Site nextSite, bottom;
		Edge newEdge;
		Vertex nextVertex, intersection;
		
		while ( pqIterator.hasNext() )
		{
			nextEvent = pqIterator.next();
			
			if ( nextEvent instanceof Site )
			{
				nextSite = (Site) nextEvent;
				// Find the first HalfEdge to the left of this site
				leftBoundary = sweepTable.leftBoundary(nextSite);
				// TODO check if this is necessary
				rightBoundary = leftBoundary.rightNeighbor;
				// 
				bottom = regionOnRight(leftBoundary);
				//
				newEdge = _space.bisect(bottom, nextSite);
				edges.add(newEdge);
				// Create
				bisector = new HalfEdge();
				bisector.edge = newEdge;
				bisector.leftRight = 0;
				// 
				sweepTable.insert(leftBoundary, bisector);
				// 
				intersection = _space.intersect(leftBoundary, bisector);
				if ( intersection != null )
				{
					// PQdelete() ???
					priorityQueueInsert(intersection, 
							_space.distance(intersection,  nextSite));
				}
				
				bisector = new HalfEdge();
				bisector.edge = newEdge;
				bisector.leftRight = 1;
				
				intersection = _space.intersect(bisector, rightBoundary);
				if ( intersection != null )
					priorityQueueInsert(intersection, 
							_space.distance(intersection,  nextSite));
				
			}
			else if ( nextEvent instanceof Vertex)
			{
				nextVertex = (Vertex) nextEvent;
				
				
				
			}
			
			
		} // End of while ( pqIterator.hasNext() )
		
		// Deal with any Edges that cross boundaries
		// TODO just those in sweepTable?
		
	}
	
	/**
	 * 
	 * @param vertex
	 * @param starOffset
	 */
	private void priorityQueueInsert(Vertex vertex, Double starOffset)
	{
		int currentIndex = pqIterator.nextIndex() - 1;
		ListIterator<ContinuousVector> insertIterator = 
									priorityQueue.listIterator(currentIndex);
		
		// TODO change .y 
		vertex.starValue = vertex.y + starOffset; 
		
		ContinuousVector next = insertIterator.next();
		
		while (insertIterator.hasNext() && _space.compare(vertex, next) > 0)
			next = insertIterator.next();
		
		insertIterator.add(vertex);
	}
	
	private Site bottomSite()
	{
		return (Site) priorityQueue.peekFirst();
	}
	
	private Site regionOnLeft(HalfEdge he)
	{
		if ( he.edge == null )
			return bottomSite();
		//TODO
		return he.edge.region[0];
	}
	
	private Site regionOnRight(HalfEdge he)
	{
		if ( he.edge == null )
			return bottomSite();
		//TODO
				return he.edge.region[0];
	}
	
	private void setEndPoint(Edge edge, int leftRight, Vertex vertex)
	{
		edge.endPoint[leftRight] = vertex;
		if ( edge.endPoint[1-leftRight] == null )
			return;
		clip(edge);
	}
	
	private void clip(Edge edge)
	{
		ContinuousVector diff = new ContinuousVector();
		LinkedList<ContinuousVector> intersections;
		diff.sendDiff(edge.endPoint[1], edge.endPoint[0]);
		for ( IsShape boundary : boundaries )
		{
			intersections = boundary.getIntersections(edge.endPoint[0], diff);
			if ( intersections == null )
				continue;
			
			Boolean outside = boundary.isOutside(edge.endPoint[0]);
			for ( ContinuousVector intersection : intersections )
			{
				
			}
			// TODO
		}
	}
}
