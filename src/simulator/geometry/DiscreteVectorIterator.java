package simulator.geometry;

public class DiscreteVectorIterator extends DiscreteVector
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public int iMin, iMax;
	
	public int jMin, jMax;
	
	public int kMin, kMax;
	
	public DiscreteVectorIterator(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
	{
		this.iMin = iMin;
		this.iMax = iMax;
		this.jMin = jMin;
		this.jMax = jMax;
		this.kMin = kMin;
		this.kMax = kMax;
		
		reset();
	}
	
	@Override
	public void reset()
	{
		set(iMin, jMin, kMin);
	}
	
	public boolean hasNext()
	{
		return ( i < iMax ) || ( j < jMax ) || ( k < kMax );
	}
	
	public boolean setNext()
	{
		if ( k < kMax )
		{
			k++;
			return true;
		}
		else if ( j < jMax )
		{
			j++;
			k = kMin;
			return true;
		}
		else if ( i < iMax )
		{
			i++;
			j = jMin;
			k = kMin;
			return true;
		}
		return false;
	}
	
	public boolean hasPrevious()
	{
		return ( i > iMin ) || ( j > jMin ) || ( k > kMin );
	}
	
	public boolean setPrevious()
	{
		if ( k > kMin )
		{
			k--;
			return true;
		}
		else if ( j > jMin )
		{
			j--;
			k = kMax;
			return true;
		}
		else if ( i > iMin )
		{
			i--;
			j = jMax;
			k = kMax;
			return true;
		}
		return false;
	}
}
