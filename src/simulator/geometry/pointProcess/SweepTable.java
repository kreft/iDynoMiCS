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

		this.minValue = shape.getMinPrimary();
		this.deltaValue = shape.getMaxPrimary() - minValue;
		
		this.hash = new HalfEdge[this.size];
		
		this.leftEnd = new HalfEdge();
		this.rightEnd = new HalfEdge();
		this.leftEnd.rightNeighbor = this.rightEnd;
		this.rightEnd.leftNeighbor = this.leftEnd;

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
		newRight.leftNeighbor = oldLeft;
		newRight.rightNeighbor = oldLeft.rightNeighbor;
		(newRight.rightNeighbor).leftNeighbor = newRight;
		oldLeft.rightNeighbor = newRight;
	}
	
	/**
	 * \brief
	 * 
	 * @param he HalfEdge to be deleted.
	 */
	public void delete(HalfEdge he)
	{
		(he.leftNeighbor).rightNeighbor = he.rightNeighbor;
		(he.rightNeighbor).leftNeighbor = he.leftNeighbor;
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
		/* Use hash table to get close to desired halfedge */
		Double temp = size*(getValue(point)-minValue)/deltaValue;
		int bucket = temp.intValue();
		/* Ensure bucket is in the range (0, this.size - 1) */
		bucket = Math.max(Math.min(bucket, size - 1), 0);
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
			do { out = out.rightNeighbor; }
			while ( out != rightEnd && isHErightOfPoint(out, point) );
			/* This is the HalfEdge immediately the right of the point, so go
			 * left one HE.
			 */
			out = out.leftNeighbor;
		}
		else
		{
			// Need a do-while in case: out == rightEnd && isLeftOfSite()
			do { out = out.leftNeighbor; }
			while ( (out != leftEnd) && isHEleftOfPoint(out, point) );
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
	private Boolean isHEleftOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		return ! isHErightOfPoint(halfEdge, point);
	}
	
	private Boolean isHErightOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		Double[] p = _space.convertToLocal(point);
		Double[] r = _space.convertToLocal(halfEdge.getRightRegion());
		Double primaryDiff = p[primary] - r[primary];
		if ( primaryDiff > 0 && halfEdge.isOnLeft() )
			return true;
		if ( primaryDiff < 0 && halfEdge.isOnRight() )
			return false;
		
		Double secondaryDiff = p[secondary] - r[secondary];
		Boolean above = true;
		// Temporary variables.
		Double a = halfEdge.edge.coefficient[0];
		Double b = halfEdge.edge.coefficient[1];
		Double c = halfEdge.edge.coefficient[2];
		Double t1, t2, t3;
		if ( halfEdge.isNearVertical() )
		{
			Boolean fast = false;
			if ( primaryDiff < 0 && b < 0.0 || primaryDiff > 0 && b >= 0.0 )
			{
				above = ( secondaryDiff >= b * primaryDiff);
				fast = above;
			}
			else
			{
				above = p[primary] + p[secondary]*b > c;
				if ( b < 0.0 )
					above = ! above;
				if ( ! above )
					fast = true;
			}
			if ( ! fast )
			{
				t1 = ExtraMath.sq(primaryDiff) - ExtraMath.sq(secondaryDiff);
				t1 *= b;
				t2 = r[primary] - _space.getPrimary(halfEdge.getLeftRegion());
				t2 *= secondaryDiff * (1.0 + ExtraMath.sq(b));
				t2 += 2.0 * primaryDiff * secondaryDiff;
				above = t1 < t2;
				if ( b < 0.0 )
					above = ! above;
			}
		}
		else
		{
			t1 = c - a*p[primary];
			t2 = p[secondary] - t1;
			t3 = t1 - r[secondary];
			above = ExtraMath.sq(t2) >
								ExtraMath.sq(primaryDiff) + ExtraMath.sq(t3); 
		}
		return halfEdge.isOnLeft() == above;
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
}