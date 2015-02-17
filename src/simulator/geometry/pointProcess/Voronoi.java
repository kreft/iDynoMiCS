package simulator.geometry.pointProcess;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.ListIterator;

import simulator.agent.LocatedGroup;
import simulator.geometry.ContinuousVector;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.shape.IsShape;
import utils.LogFile;

public class Voronoi
{
	private LinkedList<ContinuousVector> priorityQueue;
	
	private ListIterator<ContinuousVector> pqIterator;
	
	private SweepTable sweepTable;
	
	private static IsShape _space;
	
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
			//System.out.println("Added to Priority Queue: site at "+site.toString());
		}
		Collections.sort(priorityQueue, new PositionComparator());
		System.out.println("Priority Queue sorted");
		pqIterator = priorityQueue.listIterator();
		
		// Initialise the Sweep Table
		sweepTable = new SweepTable(_space, sites.size());
		
		// Initialise the list of edges.
		edges = new LinkedList<Edge>();
		
		HalfEdge leftBoundary, rightBoundary, bisector;
		Site nextSite, bottom;
		Edge newEdge;
		Vertex nextVertex, intersection;
		
		// Skip the first site, as this is the very bottom one.
		pqIterator.next();
		System.out.println("Bottom site at "+bottomSite().toString());
		
		while ( pqIterator.hasNext() )
		{
			nextEvent = pqIterator.next();
			
			if ( nextEvent instanceof Site )
			{
				System.out.println("nextEvent is a Site at "+nextEvent.toString());
				nextSite = (Site) nextEvent;
				// Find the first HalfEdge to the left of this site
				leftBoundary = sweepTable.leftBoundary(nextSite);
				// TODO check if this is necessary
				rightBoundary = leftBoundary.rightNeighbor;
				// 
				bottom = regionOnRight(leftBoundary);
				//System.out.println("nextSite at "+nextSite.toString());
				//System.out.println("bottom at "+bottom.toString());
				//
				newEdge = _space.bisect(bottom, nextSite);
				//System.out.println("newEdge "+newEdge.toString());
				edges.add(newEdge);
				// Create
				bisector = new HalfEdge(newEdge, true);
				// 
				sweepTable.insert(leftBoundary, bisector);
				// 
				//System.out.println("leftBoundary "+leftBoundary.toString());
				//System.out.println("bisector "+bisector.toString());
				intersection = _space.intersect(leftBoundary, bisector);
				if ( intersection != null )
				{
					// PQdelete() ???
					priorityQueueInsert(intersection, 
							_space.distance(intersection,  nextSite));
				}
				
				bisector = new HalfEdge(newEdge, false);
				
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
	
	
	public static class PositionComparator implements 
										java.util.Comparator<ContinuousVector> 
	{
		@Override
		public int compare(ContinuousVector point1, ContinuousVector point2) 
		{
			Double[] p1 = _space.convertToLocal(point1);
			Double[] p2 = _space.convertToLocal(point2);
			Double temp = p1[_space.getSecondary()] - p2[_space.getSecondary()];
			int out = (int) Math.signum(temp);
			if ( out == 0 )
			{
				temp = p1[_space.getPrimary()] - p2[_space.getPrimary()];
				out = (int) Math.signum(temp);
			}
			return out;
		}
	}
	
	
	/**
	 * \brief Compare two points according to the local scheme.
	 * 
	 * PriorityQueue uses the other axis to SweepTable.
	 * Primary-axis = y-axis in Fortune's paper.
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 * @see SweepTable.compare()
	 */
	public int compare(ContinuousVector point1, ContinuousVector point2) 
	{
		Double[] p1 = _space.convertToLocal(point1);
		Double[] p2 = _space.convertToLocal(point2);
		Double temp = p1[_space.getSecondary()] - p2[_space.getSecondary()];
		int out = (int) Math.signum(temp);
		if ( out == 0 )
		{
			temp = p1[_space.getPrimary()] - p2[_space.getPrimary()];
			out = (int) Math.signum(temp);
		}
		return out;
	}
	
	private Double getValue(ContinuousVector point)
	{
		return _space.getSecondary(point);
	}
	
	/**
	 * 
	 * @param vertex
	 * @param starOffset
	 */
	private void priorityQueueInsert(Vertex vertex, Double starOffset)
	{
		int currentIndex = Math.max(pqIterator.nextIndex() - 1, 0);
		ListIterator<ContinuousVector> insertIterator = 
									priorityQueue.listIterator(currentIndex);
		vertex.starValue = getValue(vertex) + starOffset;
		ContinuousVector next = insertIterator.next();
		while ( insertIterator.hasNext() && compare(vertex, next) > 0)
		{
			next = insertIterator.next();
		}
		insertIterator.add(vertex);
	}
	
	/**
	 * \brief Gets the site at the start of the PriorityQueue, i.e. the Site
	 * with the lowest PositionComparator value.
	 * 
	 * Does not set anything.
	 * 
	 * @return Site at the bottom of _space.
	 */
	private Site bottomSite()
	{
		return (Site) priorityQueue.peekFirst();
	}
	
	/**
	 * 
	 * 
	 * @param he
	 * @return
	 */
	private Site regionOnLeft(HalfEdge he)
	{
		if ( he.edge == null )
			return bottomSite();
		// Note this is the opposite outcome to regionOnRight().
		return he.isOnLeft() ?  he.getLeftRegion() : he.getRightRegion();
	}
	
	/**
	 * \brief Find the region (i.e. Site) on the right of the given HalfEdge.
	 * 
	 * Does not set anything.
	 * 
	 * @param he HalfEdge to be considered.
	 * @return Site on the right of this HalfEdge.
	 */
	private Site regionOnRight(HalfEdge he)
	{
		if ( he.edge == null )
			return bottomSite();
		// Note this is the opposite outcome to regionOnLeft().
		return he.isOnLeft() ?  he.getRightRegion() : he.getLeftRegion();
	}
	
	/**
	 * 
	 * @param edge
	 * @param leftRight
	 * @param vertex
	 */
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
