package simulator.geometry.pointProcess;

import simulator.geometry.ContinuousVector;
import simulator.geometry.shape.IsShape;
import utils.ExtraMath;

/**
 * 
 * 
 *
 */
public class SweepTable
{
	/**
	 * 
	 */
	private IsShape _space;
	
	/**
	 * 
	 */
	int size;
	
	/**
	 * 
	 */
	HalfEdge[] hash;
	
	HalfEdge leftEnd, rightEnd;
	
	/**
	 * 
	 */
	Double minValue;
	
	/**
	 * 
	 */
	Double deltaValue;
	
	private static int primary, secondary;
	
	public SweepTable(IsShape shape, int numberOfSites)
	{
		this._space = shape;
		
		primary = _space.getPrimary();
		secondary = _space.getSecondary();
		
		Double temp = 2.0 * Math.sqrt(numberOfSites + 4);
		this.size = temp.intValue();
		System.out.println("Space min "+shape.getMinPrimary()+", max "+shape.getMaxPrimary());
		this.minValue = shape.getMinPrimary();
		this.deltaValue = shape.getMaxPrimary() - minValue;
		
		this.hash = new HalfEdge[this.size];
		
		this.leftEnd = new HalfEdge();
		this.rightEnd = new HalfEdge();
		this.leftEnd.nextNeighbor = this.rightEnd;
		this.rightEnd.previousNeighbor = this.leftEnd;

		this.hash[0] = leftEnd;
		this.hash[this.size - 1] = rightEnd;
	}
	
	/**
	 * 
	 * 
	 * @param b
	 * @return
	 */
	private HalfEdge get(int b)
	{
		// If we're outside the appropriate range, return null.
		if ( b < 0 || b >= size)
			return null;
		// Find the HalfEdge corresponding to the given integer.
		HalfEdge out = hash[b];
		// If this is marked for deletion, delete it and return null.
		if ( out != null )
			if ( out.deleted )
				return (hash[b] = null);
		// Otherwise, return the HalfEdge (even if it is null).
		return out;
	}
	
	/**
	 * 
	 * 
	 * @param oldLeft
	 * @param newRight
	 */
	public void insert(HalfEdge oldLeft, HalfEdge newRight)
	{
		newRight.previousNeighbor = oldLeft;
		newRight.nextNeighbor = oldLeft.nextNeighbor;
		(newRight.nextNeighbor).previousNeighbor = newRight;
		oldLeft.nextNeighbor = newRight;
	}
	
	/**
	 * \brief
	 * 
	 * @param he HalfEdge to be deleted.
	 */
	public void delete(HalfEdge he)
	{
		(he.previousNeighbor).nextNeighbor = he.nextNeighbor;
		(he.nextNeighbor).previousNeighbor = he.previousNeighbor;
		he.deleted = true;
	}
	
	/**
	 * \brief find the HalfEdge immediately to the left of the given point.
	 * 
	 * @param point 
	 * @return
	 */
	public HalfEdge leftBoundary(ContinuousVector point)
	{
		/*
		 * Use hash table to get close to desired halfedge
		 */
		//System.out.println("size "+size+", value "+getValue(point)+", min "+minValue+", delta "+deltaValue);
		Double temp = size*(getValue(point)-minValue)/deltaValue;
		int bucket = temp.intValue();
		/* Ensure bucket is in the range (0, this.size - 1) */
		bucket = Math.max(Math.min(bucket, size - 1), 0);
		System.out.println("bucket "+bucket);
		HalfEdge out = get(bucket);
		
		/* Starting with bucket, search backwards and forwards in the hash map
		 * to find the first non-null entry. This is our initial guess.
		 */
		if ( out == null )
			for (int i = 1; i < size; i++)
			{
				if ( (out = get(bucket - i)) != null )
					break;
				if ( (out = get(bucket + i)) != null )
					break;
			}
		
		/* Linear search through the HalfEdges:
		 * 	IF: The initial guess is to the left of the Site, so keep moving
		 * 		right until the HE to the right is right of the Site.
		 *	ELSE: The initial guess is to the right of the Site, so keep
		 *		moving left until the HE is left of the Site.
		 */ 
		if (out==leftEnd || ((out!=rightEnd) && isHErightOfPoint(out, point)))
		{
			printHalfEdge("Starting left with ", out, ", going right");
			do
			{
				out = out.nextNeighbor;
				
				printHalfEdge("\tTrying ", out, "");
				if ( out != rightEnd )
					System.out.println("\t\tIs right of point? "+isHErightOfPoint(out, point));
			}
			while ( out != rightEnd && isHEleftOfPoint(out, point) );
			/* This is the HalfEdge immediately the right of the point, so go
			 * left one HE.
			 */
			out = out.previousNeighbor;
			printHalfEdge("\tGoing left one, to ", out, "");
		}
		else
		{
			printHalfEdge("Starting right with ", out, ", going left");
			// Need a do-while in case: out == rightEnd && isLeftOfSite()
			do
			{
				out = out.previousNeighbor;
				printHalfEdge("\tTrying ", out, "");
			}
			while ( out != rightEnd && isHErightOfPoint(out, point) );
			/* Nothing more to do here - we've found the HalfEdge immediately
			 * to the left of the point.
			 */
		}
		// Update hash table.
		if (bucket > 0 && bucket < this.size - 1)
			this.hash[bucket] = out;
		return out;
	}
	
	/**
	 * TODO Check!
	 * 
	 * \brief 
	 * 
	 * @param halfEdge
	 * @param point
	 * @return
	 */
	private Boolean isHErightOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		return ! isHEleftOfPoint(halfEdge, point);
	}
	
	/**
	 * 
	 * @param halfEdge
	 * @param point
	 * @return
	 */
	private Boolean isHEleftOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		/*
		 * Convert the point given, and the Site on the right of the HalfEdge,
		 * to local coordinates.
		 */
		Double[] p = _space.convertToLocal(point);
		Double[] r = _space.convertToLocal(halfEdge.getSiteAbove());
		Double primaryDiff = p[primary] - r[primary];
		/*
		 * TODO Delete above?
		 */
		Boolean rightOfSiteAbove = 
					_space.comparePrimary(point, halfEdge.getSiteAbove()) > 0;
		Boolean leftOfSiteAbove = ! rightOfSiteAbove;
		/*
		 * If the point's primary coordinate is greater than the site's, and
		 * the HE is on the left side of the Edge, then the HE is on the left
		 * of the point.
		 */
		if ( rightOfSiteAbove && halfEdge.isOutbound() )
			return true;
		/*
		 * If the point's primary coordinate is smaller than the site's, and
		 * the HE is on the right side of the Edge, then the HE is on the
		 * right of the point.
		 */
		if ( leftOfSiteAbove && halfEdge.isInbound() )
			return false;
		/*
		 * Those were the simple cases: now solve the more complicated cases.
		 */
		Double secondaryDiff = p[secondary] - r[secondary];
		Boolean pointAboveEdge = true;
		/*
		 * Unpack the edge equation to make code more readable. Names are as
		 * for edges in the plane.
		 */
		Double kv = halfEdge.edge.coefficient[1];
		/*
		 * Temporary variables to make code more readable.
		 */
		ContinuousVector pointOnEdge =
					_space.getEdgePointFromPrimary(halfEdge.edge, p[primary]);
		Double t1, t2, t3;
		if ( halfEdge.edge.isNearVertical() )
		{
			/*
			 * The edge is more parallel to the secondary axis than it is to
			 * the primary axis.
			 * 
			 * In the plane, u + kv*v = K and so v = (-1/kv)*(u - K).
			 *     If kv is negative, the slope of the edge is positive and >1
			 * (-1 < kv < 0). In other words it's pointing between 12:00 and
			 * 1:30 on a clock, or North to NW on a compass.
			 *     If kv is non-negative, the slope of the edge is negative
			 * and >1 or the line is vertical (0 <= kv < 1). In other words
			 * it's pointing between 10:30 and 12:00 on a clock, or NE to
			 * North on a compass.
			 */
			Boolean fast = false;
			if ( leftOfSiteAbove && kv < 0.0 || rightOfSiteAbove && kv >= 0.0 )
			{
				/*
				 * 
				 */
				pointAboveEdge = ( secondaryDiff >= kv * primaryDiff);
				fast = pointAboveEdge;
			}
			else
			{
				/*
				 * 
				 */
				pointAboveEdge = _space.compareSecondary(point, pointOnEdge) > 0;
				if ( kv > 0.0 )
					pointAboveEdge = ! pointAboveEdge;
				if ( ! pointAboveEdge )
					fast = true;
			}
			/*
			 * If fast is false the problem is still not solved, and so we
			 * test further.
			 */
			if ( ! fast )
			{
				t1 = ExtraMath.sq(primaryDiff) - ExtraMath.sq(secondaryDiff);
				t1 *= kv;
				t2 = r[primary] - _space.getPrimary(halfEdge.getSiteBelow());
				t2 *= secondaryDiff * (1.0 + ExtraMath.sq(kv));
				t2 += 2.0 * primaryDiff * secondaryDiff;
				pointAboveEdge = t1 < t2;
				if ( kv < 0.0 )
					pointAboveEdge = ! pointAboveEdge;
			}
		}
		else
		{
			/*
			 * TODO Explain why this is the way it is!
			 */
			pointAboveEdge = _space.distance(point, pointOnEdge) > 
					_space.distance(halfEdge.getSiteAbove(), pointOnEdge);
		}
		return halfEdge.isOutbound() == pointAboveEdge;
	}
	
	private Double getValue(ContinuousVector point)
	{
		return _space.getPrimary(point);
	}
	
	/**
	 * \brief Compare two points 
	 * 
	 * SweepTable uses the other axis to PriorityQueue.
	 * Primary-axis = x-axis in Fortune's paper.
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 * @see Voronoi.compare()
	 */
	private int compare(ContinuousVector point1, ContinuousVector point2)
	{
		Double[] p1 = _space.convertToLocal(point1);
		Double[] p2 = _space.convertToLocal(point2);
		Double temp = p1[_space.getPrimary()] - p2[_space.getPrimary()];
		int out = (int) Math.signum(temp);
		if ( out == 0 )
		{
			temp = p1[_space.getSecondary()] - p2[_space.getSecondary()];
			out = (int) Math.signum(temp);
		}
		return out;
	}
	
	public void clipAll()
	{
		// TODO
		for (HalfEdge he = leftEnd; he != rightEnd; he = he.nextNeighbor)
			he = null;
		
	}
	
	public void printSweepTable()
	{
		System.out.println("Sweep Table:");
		for ( HalfEdge he = this.leftEnd; he != null; he = he.nextNeighbor)
			System.out.println(he.toString());
	}
	
	public void printHalfEdge(String preMsg, HalfEdge he, String postMsg)
	{
		String msg = preMsg;
		if ( he.equals(this.leftEnd) )
			msg += "SweepTable leftEnd";
		else if ( he.equals(this.rightEnd) )
			msg += "SweepTable rightEnd";
		else
			msg += he.toString();
		System.out.println(msg+postMsg);
	}
}