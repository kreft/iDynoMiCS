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
import utils.ResultFile;

public class Voronoi
{
	private LinkedList<ContinuousVector> priorityQueue;
	
	private ListIterator<ContinuousVector> pqIterator;
	
	private SweepTable sweepTable;
	
	private static IsShape _space;
	
	private IsShape[] boundaries;
	
	private LinkedList<Site> sites;
	
	private LinkedList<Edge> edges;
	
	private ContinuousVector nextEvent;
	
	public Voronoi(IsShape space, LinkedList<Site> sites)
	{
		_space = space;
		/*
		 * Initialise the Priority Queue. We keep an ordered list of sites for
		 * output, to make checking/debugging easier!
		 */
		this.sites= new LinkedList<Site>();
		priorityQueue = new LinkedList<ContinuousVector>();
		for (Site site : sites)
			this.sites.add(site);
		Collections.sort(this.sites, new PositionComparator());
		for ( Site site : this.sites )
			priorityQueue.add(site);
		pqIterator = priorityQueue.listIterator();
		/*
		 * Initialise the Sweep Table.
		 */
		sweepTable = new SweepTable(_space, sites.size());
		/*
		 * Initialise the list of edges.
		 */
		edges = new LinkedList<Edge>();
		
		HalfEdge leftBoundary, rightBoundary, newHE;
		
		// TODO is this necessary?
		HalfEdge leftLeftBoundary, rightRightBoundary;
		
		Boolean isLeftHanded;
		Site nextSite, bottom, top, temp;
		Edge bisector;
		Vertex nextVertex, intersection;
		Double offset;
		/*
		 * Skip the first site, as this is the very bottom one.
		 */
		pqIterator.next();
		System.out.println("Bottom site at "+bottomSite().toString());
		
		while ( pqIterator.hasNext() )
		{
			System.out.println("");
			printPriorityQueue(pqIterator.nextIndex());
			System.out.println("");
			sweepTable.printSweepTable();
			System.out.println("");
			//printReport();
			//System.out.println("");
			nextEvent = pqIterator.next();
			
			if ( nextEvent instanceof Site )
			{
				System.out.println("nextEvent is a Site at "+nextEvent.toString());
				nextSite = (Site) nextEvent;
				/*
				 * Find the first HalfEdge to the left of this site.
				 */
				leftBoundary = sweepTable.leftBoundary(nextSite);
				/*
				 * Store the HalfEdge immediately to the right, as the left
				 * HE will change before we need it.
				 */
				rightBoundary = leftBoundary.nextNeighbor;
				/*
				 * 
				 */
				bottom = leftBoundary.getSiteAhead();
				if ( bottom == null )
					bottom = bottomSite();
				System.out.println("bottom at "+bottom.toString());
				/*
				 * Find the bisector of these two sites. Note that bottom will
				 * become bisector.region[0] and nextSite will become 
				 * bisector.region[1].
				 */
				bisector = _space.bisect(bottom, nextSite);
				edges.add(bisector);
				/*
				 * Create a HalfEdge on the left of this bisector.
				 */
				newHE = new HalfEdge(bisector, true);
				/*
				 * 
				 */
				System.out.println("Inserting");
				System.out.println("\tbisector: "+newHE.toString());
				System.out.println("to the right of");
				System.out.println("\tleftBoundary: "+leftBoundary.toString());
				sweepTable.insert(leftBoundary, newHE);
				intersection = intersect(leftBoundary, newHE);
				if ( intersection != null )
				{
					// PQdelete() ???
					System.out.println("Vertex found at "+intersection.toString()+" [A]");
					offset = _space.distance(intersection,  nextSite);
					priorityQueueInsert(intersection, offset);
				}
				/*
				 * Now look at the right-hand HalfEdge of newEdge.
				 */
				leftBoundary = newHE;
				newHE = new HalfEdge(bisector, false);
				System.out.println("Inserting");
				System.out.println("\tbisector: "+newHE.toString());
				System.out.println("to the right of");
				System.out.println("\tleftBoundary: "+leftBoundary.toString());
				sweepTable.insert(leftBoundary, newHE);
				/*
				 * See if this intersects rightBoundary.
				 */
				intersection = intersect(newHE, rightBoundary);
				if ( intersection != null )
				{
					System.out.println("Vertex found at "+intersection.toString()+" [B]");
					offset = _space.distance(intersection,  nextSite);
					priorityQueueInsert(intersection, offset);
				}
			}
			else if ( nextEvent instanceof Vertex)
			{
				System.out.println("nextEvent is a Vertex at "+nextEvent.toString());
				nextVertex = (Vertex) nextEvent;
				/*
				 * Find the first HalfEdge to the left of this vertex.
				 */
				leftBoundary = sweepTable.leftBoundary(nextVertex);
				/*
				 * 
				 */
				leftLeftBoundary = leftBoundary.previousNeighbor;
				/*
				 * Store the HalfEdge immediately to the right, as the left
				 * HE will change before we need it.
				 */
				rightBoundary = leftBoundary.nextNeighbor;
				/*
				 * 
				 */
				rightRightBoundary = rightBoundary.nextNeighbor;
				/*
				 * Get the site to the left of the HalfEdge on the left.
				 */
				bottom = leftBoundary.getSiteBehind();
				if ( bottom == null )
					bottom = bottomSite();
				System.out.println("bottom at "+bottom.toString());
				/*
				 * Get the site to the right of the HalfEdge on the right.
				 */
				top = rightBoundary.getSiteAhead();
				if ( top == null )
					top = bottomSite();
				// TODO set vertex number?
				/* 
				 * 
				 */
				System.out.println("HalfEdge on left: "+leftBoundary.toString());
				setEndPoint(leftBoundary.edge, leftBoundary.isOutbound(), nextVertex);
				System.out.println("HalfEdge on right: "+rightBoundary.toString());
				setEndPoint(rightBoundary.edge, rightBoundary.isOutbound(), nextVertex);
				/*
				 * Assume the new bisector will be left-handed.
				 */
				isLeftHanded = true;
				/*
				 * If bottom is higher than top, switch them and change the
				 * handedness of the new bisector to right.
				 */
				if ( _space.compareSecondary(bottom, top) > 0 )
				{
					temp = bottom;
					bottom = top;
					top = temp;
					isLeftHanded = false;
				}
				bisector = _space.bisect(bottom, top);
				edges.add(bisector);
				newHE = new HalfEdge(bisector, isLeftHanded);
				sweepTable.insert(leftBoundary, newHE);
				System.out.println("newEdge "+bisector.toString());
				setEndPoint(bisector, ! isLeftHanded, nextVertex);
				
				intersection = intersect(leftLeftBoundary, newHE);
				if ( intersection != null )
				{
					System.out.println("Vertex found at "+intersection.toString()+" [C]");
					offset = _space.distance(intersection,  bottom);
					priorityQueueInsert(intersection, offset);
				}
				/*
				 * TODO Check why it's the distance to the bottom that's used
				 * here - the top seems to make more sense.
				 */
				intersection = intersect(newHE, rightRightBoundary);
				if ( intersection != null )
				{
					System.out.println("Vertex found at "+intersection.toString()+" [D]");
					offset = _space.distance(intersection,  bottom);
					priorityQueueInsert(intersection, offset);
				}
			}
			
			
		} // End of while ( pqIterator.hasNext() )
		
		
		
		// Deal with any Edges that cross boundaries
		// TODO just those in sweepTable?
		//for ( Edge edge : edges )
		//	clip(edge);
		
		System.out.println("");
		System.out.println("------------------- FINISHED! -------------------");
		System.out.println("");
		printPriorityQueue(-1);
		System.out.println("");
		sweepTable.printSweepTable();
		System.out.println("");
		printReport();
	}
	
	/**
	 * \brief Comparator for points, so they can be lexicographically ordered.
	 * 
	 * Fortune's paper states that, for the purpose of the priority queue,
	 * points in the xy-plane should first be ordered be y and then by x. I.e.
	 * p < q if py < qy, or if py = qy and px < px. Since here we are less
	 * restricted in the surfaces we work over, we order first according to
	 * the "secondary" value of points (equivalent to their py, etc) and then
	 * according to their "primary" value (equivalent to px, etc).
	 */
	public static class PositionComparator implements 
										java.util.Comparator<ContinuousVector> 
	{
		@Override
		public int compare(ContinuousVector point1, ContinuousVector point2) 
		{
			return _space.compareSecondary(point1, point2);
		}
	}
	
	
	private Double getValue(ContinuousVector point)
	{
		return _space.getSecondary(point);
	}
	
	private Vertex intersect(HalfEdge he1, HalfEdge he2)
	{
		
		Vertex out = _space.intersect(he1.edge, he2.edge);
		if ( out == null )
			return out;
		/*
		 * Find the HalfEdge with the "greater" region on the right.
		 * 
		 * TODO Check!!!
		 */
		HalfEdge temp;
		if ( _space.compareSecondary(he1.getSiteAbove(), he2.getSiteAbove()) > 0 )
			temp = he1;
		else
			temp = he2;
		/*
		 * 
		 */
		Boolean rightOfSite = getValue(out) >= getValue(temp.getSiteAbove());
		/*
		 * 
		 */
		if ( rightOfSite == temp.isOutbound() )
			return null;
		/*
		 * Once we've passed all these checks, the vertex can be returned.
		 */
		return out;
	}
	
	/**
	 * 
	 * TODO Tidy up
	 * 
	 * @param vertex
	 * @param starOffset
	 */
	private void priorityQueueInsert(Vertex vertex, Double starOffset)
	{
		int currentIndex = pqIterator.previousIndex();
		//System.out.println("index of event = "+currentIndex);
		ContinuousVector current = pqIterator.previous();
		ContinuousVector next;
		vertex.starValue = getValue(vertex) + starOffset;
		/*
		 * Search the priority queue for a position with
		 * _space.compareSecondary(vertex, next) < 0
		 * 
		 * If _space.compareSecondary(vertex, current) < 0
		 * we take the next position
		 */
		do { next = pqIterator.next(); }
		while ( pqIterator.hasNext() && _space.compareSecondary(vertex, next) > 0 );
		/*
		 * Once we've found the right place, add the vertex.
		 */
		pqIterator.add(vertex);
		//currentIndex = pqIterator.previousIndex();
		//System.out.println("index of inserted vertex = "+currentIndex);
		/*
		 * Get back to where we were before.
		 */
		//while ( (! next.equals(current)) && pqIterator.hasPrevious() )
		//	next = pqIterator.previous();
		pqIterator = priorityQueue.listIterator(currentIndex+1);
		
	}
	
	/**
	 * Useful for checking/debugging.
	 * 
	 * @param index Position in the priority queue that you want to flag.
	 */
	@SuppressWarnings("unused")
	private void printPriorityQueue(int index)
	{
		String msg;
		ContinuousVector p;
		for ( int i = 0; i < priorityQueue.size(); i++ )
		{
			msg = "\t";
			if ( i == index )
				msg += "-> ";
			else
				msg += "   ";
			p = priorityQueue.get(i);
			if ( p instanceof Site )
				msg += "S ";
			else if ( p instanceof Vertex )
				msg += "V ";
			else
				msg += "p ";
			System.out.println(msg+p.toString());
		}
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
	 * TODO Check inner-outer
	 * 
	 * 
	 * 
	 * @param edge
	 * @param isOutbound
	 * @param vertex
	 */
	private void setEndPoint(Edge edge, Boolean isOutbound, Vertex vertex)
	{
		String msg = "Setting ";
		msg += isOutbound ? "inner" : "outer";
		msg += " vertex of edge to ";
		System.out.println(msg+vertex.toString());
		edge.setEndPoint(vertex, isOutbound);
		if ( edge.areEndPointsSet() )
			clip(edge);
	}
	
	/**
	 * \brief Clip an Edge that has both endPoints set.
	 * 
	 * Takes care of endPoints that go outside the boundaries.
	 * 
	 * @param edge
	 */
	private void clip(Edge edge)
	{
		ContinuousVector diff = new ContinuousVector();
		LinkedList<ContinuousVector> intersections;
		if ( ! edge.areEndPointsSet() )
		{
			//ContinuousVector[] temp = _space.
			
		}
		diff.sendDiff(edge.endPoint[1], edge.endPoint[0]);
		for ( IsShape boundary : boundaries )
		{
			intersections = boundary.getIntersections(edge.endPoint[0], diff);
			if ( intersections == null )
				continue;
			/*
			 * TODO Currently only handles one intersection, should make more
			 * robust.
			 */
			if ( boundary.isOutside(edge.endPoint[0]) )
				edge.endPoint[0].set(intersections.getFirst());
			else
				edge.endPoint[1].set(intersections.getFirst());
		}
	}
	
	/**
	 * Currently just prints the report to screen.
	 * 
	 * TODO make this write the report to a ResultFile
	 */
	public void printReport()
	{
		/*
		 * This will be our general-purpose buffer.
		 */
		StringBuffer textBuffer = new StringBuffer();
		/*
		 * First, fill the buffer with information about the space.
		 */
		_space.writeShapeInformation(textBuffer);
		/*
		 * Now create a list of the sites.
		 */
		textBuffer.append("<sites header=\""+_space.getSitesHeader()+"\">\n");
		for ( Site site : sites )
			textBuffer.append(_space.getVectorOutput(site)+";\n");
		textBuffer.append("<sites/>\n");
		/*
		 * Next, the edges.
		 */
		/*
		textBuffer.append("<edges header=\""+_space.getEdgesHeader()+"\">\n");
		for ( Edge edge : edges )
		{
			if ( edge.endPoint[0] == null || edge.endPoint[1] == null)
				continue;
			textBuffer.append(_space.getVectorOutput(edge.endPoint[0])+","+
							  _space.getVectorOutput(edge.endPoint[1])+";\n");
		}
		textBuffer.append("<edges/>\n");
		*/
		textBuffer.append("<edges header=\"ku,kv,K,u0,v0,u1,v1\">\n");
		Double[] p;
		for ( Edge edge : edges )
		{
			textBuffer.append(edge.coefficient[0]+","+
								edge.coefficient[1]+","+
								edge.coefficient[2]+",");
			if (edge.endPoint[0] == null)
				textBuffer.append("null,null,");
			else
			{
				p = _space.convertToLocal(edge.endPoint[0]);
				textBuffer.append(p[0]+","+p[1]+",");
			}
			if (edge.endPoint[1] == null)
				textBuffer.append("null,null");
			else
			{
				p = _space.convertToLocal(edge.endPoint[1]);
				textBuffer.append(p[0]+","+p[1]);
			}
			textBuffer.append(";\n");
		}
		textBuffer.append("<edges/>\n");
		System.out.print(textBuffer.toString());
	}
}
