package simulator.geometry.pointProcess;

import java.util.Collections;
import java.util.LinkedList;
import java.util.ListIterator;

import simulator.geometry.ContinuousVector;
import simulator.geometry.pointProcess.HalfEdge;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.shape.IsShape;
import simulator.geometry.shape.Planar;
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
		
		Boolean isOutbound;
		Site nextSite, bottom, top, temp;
		Edge bisector;
		Vertex nextVertex, intersection;
		Double offset;
		/*
		 * Skip the first site, as this is the very bottom one.
		 */
		pqIterator.next();
		System.out.println("Bottom site at "+bottomSite().toString());
		try
		{
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
				//leftBoundary = sweepTable.leftBoundary(nextSite);
				leftBoundary = sweepTable.halfEdgeImmediatelyBehind(nextSite);
				System.out.println("leftBoundary is "+leftBoundary.toString());
				/*
				 * Store the HalfEdge immediately to the right, as the left
				 * HE will change before we need it.
				 */
				rightBoundary = leftBoundary.nextNeighbor;
				System.out.println("rightBoundary is "+rightBoundary.toString());
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
				sweepTable.insert(leftBoundary, newHE);
				System.out.println("\nTrying to intersect [A]");
				System.out.println(leftBoundary.toString());
				System.out.println(newHE.toString());
				intersection = intersect(leftBoundary, newHE);
				if ( intersection != null )
				{
					// PQdelete() ???
					System.out.println("Vertex found at "+
											intersection.toString()+" [A]");
					offset = _space.distance(intersection,  nextSite);
					priorityQueueInsert(intersection, offset);
					intersection.previousHE = leftBoundary;
				}
				else
					System.out.println("No vertex found [A]");
				/*
				 * Now look at the right-hand HalfEdge of newEdge.
				 */
				leftBoundary = newHE;
				newHE = new HalfEdge(bisector, false);
				sweepTable.insert(leftBoundary, newHE);
				/*
				 * See if this intersects rightBoundary.
				 */
				System.out.println("\nTrying to intersect [B]");
				System.out.println(newHE.toString());
				System.out.println(rightBoundary.toString());
				intersection = intersect(newHE, rightBoundary);
				if ( intersection != null )
				{
					System.out.println("Vertex found at "+
											intersection.toString()+" [B]");
					offset = _space.distance(intersection,  nextSite);
					priorityQueueInsert(intersection, offset);
					intersection.previousHE = newHE;
				}
				else
					System.out.println("No vertex found [B]");
					
			}
			else if ( nextEvent instanceof Vertex)
			{
				System.out.println("nextEvent is a Vertex at "+nextEvent.toString());
				nextVertex = (Vertex) nextEvent;
				/*
				 * Find the first HalfEdge to the left of this vertex.
				 */
				leftBoundary = nextVertex.previousHE;
				System.out.println("HE immediately to left is "+
													leftBoundary.toString());
				/*
				 * 
				 */
				leftLeftBoundary = leftBoundary.previousNeighbor;
				/*
				 * Store the HalfEdge immediately to the right, as the left
				 * HE will change before we need it.
				 */
				rightBoundary = leftBoundary.nextNeighbor;
				System.out.println("HE immediately to right is "+
													rightBoundary.toString());
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
				System.out.println("top at "+top.toString());
				// TODO set vertex number?
				/* 
				 * 
				 */
				//System.out.println("HalfEdge on left: "+leftBoundary.toString());
				setEndPoint(leftBoundary.edge, false, nextVertex);
				//System.out.println("HalfEdge on right: "+rightBoundary.toString());
				setEndPoint(rightBoundary.edge, true, nextVertex);
				/*
				 * Assume the new bisector will be out-bound.
				 */
				isOutbound = true;
				/*
				 * If bottom is higher than top, switch them and make the
				 * new half-edge in-bound.
				 */
				if ( _space.compareSecondary(bottom, top) > 0 )
				{
					temp = bottom;
					bottom = top;
					top = temp;
					isOutbound = false;
				}
				bisector = _space.bisect(bottom, top);
				edges.add(bisector);
				newHE = new HalfEdge(bisector, isOutbound);
				sweepTable.insert(leftBoundary, newHE);
				System.out.println("newEdge "+bisector.toString());
				setEndPoint(bisector, ! isOutbound, nextVertex);
				
				intersection = intersect(leftLeftBoundary, newHE);
				if ( intersection != null )
				{
					System.out.println("\nVertex found at "+
											intersection.toString()+" [C]");
					offset = _space.distance(intersection,  bottom);
					priorityQueueInsert(intersection, offset);
					intersection.previousHE = leftLeftBoundary;
				}
				/*
				 * TODO Check why it's the distance to the bottom that's used
				 * here - the top seems to make more sense.
				 */
				intersection = intersect(newHE, rightRightBoundary);
				if ( intersection != null )
				{
					System.out.println("\nVertex found at "+
											intersection.toString()+" [D]");
					offset = _space.distance(intersection,  bottom);
					priorityQueueInsert(intersection, offset);
					intersection.previousHE = newHE;
				}
			}
			
			
		} // End of while ( pqIterator.hasNext() )
		}
		catch (Exception e)
		{
			
		}
		
		
		// Deal with any Edges that cross boundaries
		// TODO just those in sweepTable?
		try
		{
		for ( Edge edge : edges )
			clip(edge);
		}
		catch (Exception e)
		{
			
		}
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
	 * \brief Insert a Vertex into the Priority Queue
	 * 
	 * @param vertex
	 * @param starOffset
	 */
	private void priorityQueueInsert(Vertex vertex, Double starOffset)
	{
		int currentIndex = pqIterator.nextIndex();
		vertex.starValue = getValue(vertex) + starOffset;
		/*
		 * Search the priority queue for a position with
		 * _space.compareSecondary(vertex, next) < 0
		 * 
		 * If _space.compareSecondary(vertex, current) < 0
		 * we take the next position
		 */
		ContinuousVector next = pqIterator.previous();
		do { next = pqIterator.next(); }
		while ( pqIterator.hasNext() && 
								_space.compareSecondary(vertex, next) > 0 );
		/*
		 * Once we've found the right place, add the vertex.
		 */
		pqIterator.add(vertex);
		/*
		 * Get back to where we were before.
		 */
		pqIterator = priorityQueue.listIterator(currentIndex);
		
	}
	
	/**
	 * Useful for checking/debugging.
	 * 
	 * @param index Position in the priority queue that you want to flag.
	 */
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
	 * \brief Set an end point of the given Edge to be the given Vertex.
	 * 
	 * Also clips the line if both end points have been set.
	 * 
	 * @param edge	The Edge to have its end point set.
	 * @param setOuter	Whether to set the outer end point (true) or the inner
	 * end point (false).
	 * @param vertex	The Vertex to be set as an end point.
	 */
	private void setEndPoint(Edge edge, Boolean setOuter, Vertex vertex)
	{
		System.out.print("setEndPoint: ");
		System.out.print(setOuter ? "outer" : "inner");
		System.out.print(" vertex "+vertex.toString()+" for edge ");
		System.out.println(edge.toString());
		
		edge.setEndPoint(vertex, setOuter);
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
		System.out.println("Voronoi.clip(Edge) looking at "+edge.toString());
		/*
		 * Build a list of all points on the shape we need to consider. This
		 * obviously includes the endPoints already set, if they exist.
		 */
		LinkedList<ContinuousVector> allToConsider =
										new LinkedList<ContinuousVector>();
		for ( Vertex ep : edge.endPoint )
			if ( ep != null )
				allToConsider.add(ep);
		/*
		 * Now add all intersections with boundaries of the shape.
		 */
		ContinuousVector diff = new ContinuousVector();
		diff.sendDiff(edge.endPoint[1], edge.endPoint[0]);
		LinkedList<ContinuousVector> intersections;
		for ( IsShape boundary : _space.getBoundaries() )
		{ 
			intersections = boundary.getIntersections(edge.endPoint[0], diff);
			if ( intersections.isEmpty() )
				continue;
			allToConsider.addAll(intersections);
		}
		/*
		 * Remove any that are outside any other boundary.
		 */
		LinkedList<ContinuousVector> toRemove = 
										new LinkedList<ContinuousVector>();
		//System.out.println("Considering...");
		for ( ContinuousVector v : allToConsider )
		{
			//System.out.println("\n\t"+v.toString());
			sLoop: for ( IsShape boundary : _space.getBoundaries() )
			{
				//System.out.println("\t\t"+boundary.toString());
				if ( boundary instanceof Planar )
				{
					Double angle = ((Planar) boundary).getNormalContinuous().
							cosAngle(((Planar) boundary).getRelativePosition(v));
					//System.out.println("\t\t\t"+angle);
				}
				if ( boundary.isOutside(v) )
				{
					//System.out.println("\tOUTSIDE!");
					toRemove.add(v);
					break sLoop;
				}
			}
		}
		allToConsider.removeAll(toRemove);
		/*
		 * Sort this list lexicographically, i.e. the same way as we did with
		 * the sites at the beginning.
		 */
		allToConsider.sort(new PositionComparator());
		switch ( allToConsider.size() )
		{
		case 2:
			edge.endPoint[0] = new Vertex();
			edge.endPoint[1] = new Vertex();
			edge.endPoint[0].set(allToConsider.getFirst());
			edge.endPoint[1].set(allToConsider.getLast());
			break;
		case 3:
			if ( allToConsider.get(1).equals(edge.endPoint[0]) )
			{
				edge.endPoint[1] = new Vertex();
				edge.endPoint[1].set(allToConsider.getLast());
			}
			else if ( allToConsider.get(1).equals(edge.endPoint[1]) )
			{
				edge.endPoint[0] = new Vertex();
				edge.endPoint[0].set(allToConsider.getFirst());
			}
			else
			{
				System.out.println("Voronoi.clip(Edge) has broken!");
				for ( ContinuousVector v : allToConsider )
					System.out.println("\t"+v.toString());
				System.exit(-1);
			}
			break;
		case 4:
			edge.endPoint[0] = new Vertex();
			edge.endPoint[1] = new Vertex();
			edge.endPoint[0].set(allToConsider.get(1));
			edge.endPoint[1].set(allToConsider.get(2));
			break;

		default:
			System.out.println("Voronoi.clip(Edge) has broken!");
			for ( ContinuousVector v : allToConsider )
				System.out.println("\t"+v.toString());
			System.exit(-1);
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
			try
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
			catch ( Exception e)
			{
				
			}
		}
		textBuffer.append("<edges/>\n");
		System.out.print(textBuffer.toString());
	}
}
