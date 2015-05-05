package utils.nDimensionalArray;



public class NDimensionalArray
{
	
	protected int[] dimensionSizes;
	
	protected int nDims;
	
	protected int totalLength;
	
	protected int[] padding;
	
	
	public NDimensionalArray(int[] dimensionSizes) throws Exception
	{
		this.dimensionSizes = dimensionSizes;
		
		this.nDims = dimensionSizes.length;
		if ( this.nDims == 0 )
			throw new Exception();
		
		this.totalLength = 1;
		padding = new int[this.nDims];
		for (int d : dimensionSizes)
		{
			if ( d <= 0 )
				throw new Exception();
			this.totalLength *= d;
			padding[d] = 0;
		}
	}
	
	public int getNDimensions()
	{
		return this.nDims;
	}
	
	public int[] getDimensionSizes()
	{
		return this.dimensionSizes;
	}
	
	/**
	 * 
	 * @param other
	 * @return
	 */
	protected Boolean areDimensionsSame(NDimensionalArray other)
	{
		/*
		 * If the other array has a different number of dimensions, then it's
		 * definitely different.
		 */
		if ( this.nDims != other.getNDimensions() )
			return false;
		/*
		 * Check for any differences between individual dimensions.
		 */
		for ( int d = 0; d < this.nDims; d++ )
			if ( this.dimensionSizes[d] != other.getDimensionSizes()[d] )
				return false;
		/*
		 * If they're all the same, then the arrays have the same dimensions.
		 */
		return true;
	}
	
	protected Boolean isInPadding(int[] position)
	{
		for ( int d = 0; d < this.nDims; d++)
		{
			if ( position[d] < this.padding[d] ||
					position[d] >= this.dimensionSizes[d] - this.padding[d])
			{
				return true;
			}
		}
		return false;
	}
	
	
	/*************************************************************************
	 * INDEXING
	 * (Could be moved to separate class for different space-filling curves)
	 ************************************************************************/
	
	/**
	 * \brief
	 * 
	 * @param position
	 * @return
	 * @throws Exception
	 */
	protected int findIndex(int[] position) throws Exception
	{
		/*
		 * Check that the position array has the correct number of dimensions.
		 */
		if ( this.nDims != position.length )
			throw new Exception();
		/*
		 * Calculate the index.
		 */
		int index = 0;
		int buffer = 1;
		for (int d = 0; d < this.nDims; d++)
		{
			index += buffer * position[d];
			buffer *= this.dimensionSizes[d];
		}
		return index;
	}
	
	/**
	 * \brief 
	 * 
	 * @param index
	 * @return
	 * @throws Exception
	 */
	protected int[] findPosition(int index) throws Exception
	{
		/*
		 * Check the index is possible.
		 */
		if ( index >= this.totalLength )
			throw new Exception();
		/*
		 * Find the position.
		 */
		int[] position = new int[this.nDims];
		int buffer = this.totalLength;
		for (int d = this.nDims - 1; d >= 0; d--)
		{
			buffer = buffer/ this.dimensionSizes[d];
			position[d] = index / buffer;
			index = index % buffer;
		}
		return position;
	}
	
	
	
	
}
