package simulator.geometry.pointProcess;

import simulator.geometry.ContinuousVector;
import simulator.geometry.shape.IsShape;

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
	IsShape space;
	
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
	
	public SweepTable(IsShape space, int numberOfSites,
										Double minValue, Double maxValue)
	{
		this.space = space;

		Double temp = 2 * Math.sqrt(numberOfSites);
		this.size = temp.intValue();

		this.minValue = minValue;
		this.deltaValue = maxValue - minValue;

		for (int i = 0; i < this.size; i++)
			this.hash[i] = null;

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
		if ( b < 0 || b >= this.size)
			return null;
		// Find the HalfEdge corresponding to the given integer.
		HalfEdge out = this.hash[b];
		// If this is marked for deletion, delete it and return null.
		if ( out != null && out.deleted )
			return (this.hash[b] = null);
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
		// TODO using p.x might not be correct!!!
		Double temp = this.size * (point.x - this.minValue)/this.deltaValue;
		int bucket = temp.intValue();
		/* Ensure bucket is in the range (0, this.size - 1) */
		bucket = Math.max(Math.min(bucket, this.size - 1), 0);
		HalfEdge out = get(bucket);
		
		/* Starting with bucket, search backwards and forwards in the hash map
		 * to find the first non-null entry. This is our initial guess.
		 */
		if ( out == null )
			for (int i = 1; i < this.size; i++)
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
		if ( out == this.leftEnd || 
					((out != this.rightEnd) && isHEleftOfPoint(out, point)))
		{
			while ( (out.rightNeighbor != this.rightEnd) && 
								isHEleftOfPoint(out.rightNeighbor, point))
				out = out.rightNeighbor;
		}
		else
		{
			// Need a do-while in case: out == rightEnd && isLeftOfSite()
			do
			{
				out = out.leftNeighbor;
			} while ( (out != leftEnd) && isHErightOfPoint(out, point) );
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
		int temp = space.compare(point, halfEdge.edge.region[1]);
		
		if ( temp > 0 && halfEdge.leftRight == 0)
			return true;
		
		if ( temp < 0 && halfEdge.leftRight == 1)
			return false;
		
		return true;
	}
	
	private Boolean isHErightOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		return ( ! isHEleftOfPoint(halfEdge, point) );
	}
}
