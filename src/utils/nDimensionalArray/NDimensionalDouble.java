package utils.nDimensionalArray;

import java.util.function.BiFunction;
import java.util.function.DoubleFunction;
import java.util.function.Function;

import utils.ExtraMath;

public class NDimensionalDouble extends NDimensionalArray
{
	private Double[] values;
	
	private Double temp;
	
	public NDimensionalDouble(int[] dimensionSizes) throws Exception
	{
		super(dimensionSizes);
		
		this.values = new Double[this.totalLength];
		this.resetAllToZero();
	}
	
	/*************************************************************************
	 * FUNCTIONAL METHODS
	 ************************************************************************/
	
	/**
	 * \brief Apply a given function to the value at the given position.
	 * 
	 * @param position
	 * @param f
	 * @throws Exception thrown if position given is outside the array.
	 */
	public void applyToValue(int[] position, DoubleFunction<Double> f)
															throws Exception
	{
		int index = findIndex(position);
		this.values[index] = (Double) f.apply((Double) this.values[index]);
	}
	
	/**
	 * \brief Apply a given function to all values.
	 * 
	 * @param f Function to apply to all values.
	 */
	public void applyToAll(DoubleFunction<Double> f, Boolean excludePadding)
	{
		try
		{
			for (int i = 0; i < this.totalLength; i++)
			{
				if ( excludePadding && isInPadding(findPosition(i)) )
					continue;
				this.values[i] = (Double) f.apply((Double) this.values[i]);
			}
		}
		catch (Exception e)
		{
			/*
			 * This shouldn't ever happen, since we are looping through the
			 * indices rather than choosing positions externally.
			 */
		}
	}
	
	/**
	 * \brief Apply a given function to all values.
	 * 
	 * @param f Function to apply to all values.
	 */
	public void getFromAll(BiFunction<Double, Double, Double> f,
														Boolean excludePadding)
	{
		try
		{
			for (int i = 0; i < this.totalLength; i++)
			{
				if ( excludePadding && isInPadding(findPosition(i)) )
					continue;
				temp = (Double) f.apply((Double) this.values[i], (Double) temp);
			}
		}
		catch (Exception e)
		{
			/*
			 * This shouldn't ever happen, since we are looping through the
			 * indices rather than choosing positions externally.
			 */
		}
	}
	
	/**
	 * \brief 
	 * 
	 * @param f
	 * @param excludePadding
	 * @return
	 * @throws Exception
	 */
	public Double getProcessedSum(DoubleFunction<Double> f,
														Boolean excludePadding)
	{
		Double out = 0.0;
		try
		{
			for (int i = 0; i < this.totalLength; i++)
			{
				if ( excludePadding && isInPadding(findPosition(i)) )
					continue;
				out += (Double) f.apply((Double) this.values[i]);
			}
		}
		catch (Exception e)
		{
			/*
			 * This shouldn't ever happen, since we are looping through the
			 * indices rather than choosing positions externally.
			 */
		}
		return out;
	}
	
	/**
	 * \brief 
	 * 
	 * @param f
	 * @param excludePadding
	 * @return
	 */
	public Boolean checkAllTrue(Function<Double, Boolean> f,
														Boolean excludePadding)
	{
		try
		{
			for (int i = 0; i < this.totalLength; i++)
			{
				if ( excludePadding && isInPadding(findPosition(i)) )
					continue;
				if ( ! f.apply(this.values[i]) )
					return false;
			}
		}
		catch (Exception e)
		{
			/*
			 * This shouldn't ever happen, since we are looping through the
			 * indices rather than choosing positions externally.
			 */
		}
		return true;
	}
	
	/**
	 * TODO Assumes both arrays are using the same indexing method
	 * 
	 * @param f
	 * @param other
	 * @param excludePadding
	 * @throws Exception
	 */
	public void applyToAll(BiFunction<Double, Double, Double> f,
			NDimensionalDouble other, Boolean excludePadding) throws Exception
	{
		/*
		 * Safety checking first.
		 */
		if ( ! areDimensionsSame(other) )
			throw new Exception();
		for (int i = 0; i < this.totalLength; i++)
		{
			if ( excludePadding && isInPadding(findPosition(i)) )
				continue;
			this.values[i] = (Double) f.apply(this.values[i], other.getValue(i));
		}
	}
	
	/**
	 * TODO Assumes both arrays are using the same indexing method
	 * 
	 * @param f
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception
	 */
	public Double getProcessedBiSum(BiFunction<Double, Double, Double> f,
			NDimensionalDouble other, Boolean excludePadding) throws Exception
	{
		/*
		 * Safety checking first.
		 */
		if ( ! areDimensionsSame(other) )
			throw new Exception();
		/*
		 * Now calculate the new sum.
		 */
		Double out = 0.0;
		for (int i = 0; i < this.totalLength; i++)
		{
			if ( excludePadding && isInPadding(findPosition(i)) )
				continue;
			out += (Double) f.apply(this.values[i], other.getValue(i));
		}
		return out;
	}
	
	/**
	 * TODO Assumes both arrays are using the same indexing method
	 * 
	 * @param f
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception
	 */
	public Boolean checkAllTrue(BiFunction<Double, Double, Boolean> f,
			NDimensionalDouble other, Boolean excludePadding) throws Exception
	{
		/*
		 * Safety checking first.
		 */
		if ( ! areDimensionsSame(other) )
			throw new Exception();
		/*
		 * Now check condition is met for all pairs of corresponding values.
		 */
		for (int i = 0; i < this.totalLength; i++)
		{
			if ( excludePadding && isInPadding(findPosition(i)) )
				continue;
			if ( ! f.apply(this.values[i], other.getValue(i)) )
				return false;
		}
		return true;
	}
	
	
	/*************************************************************************
	 * USEFUL METHODS
	 ************************************************************************/
	
	/**
	 * \brief 
	 * 
	 * @param index
	 * @return
	 */
	public Double getValue(int index)
	{
		return this.values[index];
	}
	
	/**
	 * \brief
	 * 
	 * @param position
	 * @return
	 * @throws Exception
	 */
	public Double getValue(int[] position) throws Exception
	{
		return getValue(this.findIndex(position));
	}
	
	/**
	 * \brief 
	 * 
	 * @param position
	 * @param newValue
	 * @throws Exception
	 */
	public void setValue(int[] position, Double newValue) throws Exception
	{
		this.values[this.findIndex(position)] = newValue;
	}
	
	/**
	 * \brief Sets all values to that given.
	 * 
	 * @param newValue Double value to set all values to.
	 */
	public void setAllTo(Double newValue, Boolean excludePadding)
	{
		this.applyToAll((value) -> {return newValue;}, excludePadding);
	}
	
	/**
	 * \brief Sets all values to zero.
	 * 
	 * Includes padding, if any.
	 */
	public void resetAllToZero()
	{
		this.setAllTo(0.0, true);
	}
	
	public Double max(Boolean excludePadding)
	{
		temp = Double.NEGATIVE_INFINITY;
		this.applyToAll(
				(value) -> {temp = Math.max(temp,  value); return value;},
										excludePadding);
		return temp;
	}
	
	/**
	 * \brief 
	 * 
	 * @param multiplier
	 * @param excludePadding
	 */
	public void times(Double multiplier, Boolean excludePadding)
	{
		this.applyToAll((value) -> {return multiplier * value;},
															excludePadding);
	}
	
	/**
	 * \brief If any values are negative, sets them to zero.
	 */
	public void ensureAllNonNegative()
	{
		this.applyToAll((value) -> {return Math.max(value,  0.0);}, false);
	}
	
	/**
	 * \brief 
	 * 
	 * @param excludePadding 
	 * @return Boolean stating whether array contains any negative values (False)
	 * or all values are non-negative (true).
	 */
	public Boolean areAllNonNegative(Boolean excludePadding)
	{
		try
		{
			return this.checkAllTrue((value) -> {return value >= 0.0;},
					excludePadding);
		}
		catch (Exception e)
		{
			return false;
		}
	}
	
	/**
	 * \brief
	 * 
	 * @param excludePadding
	 * @return
	 */
	public Boolean areAllZero(Boolean excludePadding)
	{
		return this.checkAllTrue((value) -> {return value == 0.0;},
															excludePadding);
	}
	
	/**
	 * \brief Checks whether the array contains any NaN or infinities.
	 * 
	 * @return
	 */
	public Boolean isFinite(Boolean excludePadding)
	{
		return this.checkAllTrue((value) -> {return Double.isFinite(value);},
															excludePadding);
	}
	
	
	public Double norm(Boolean excludePadding)
	{
		Double sumSq = this.getProcessedSum(
					(value) -> {return ExtraMath.sq(value);}, excludePadding);
		return Math.sqrt(sumSq);
	}
	
	/**
	 * \brief Sets the norm to 1.0
	 * 
	 * @param newNorm
	 * @param excludePadding
	 */
	public void normalise(Double newNorm, Boolean excludePadding)
	{
		Double oldNorm = this.norm(excludePadding);
		if ( oldNorm != 0.0 || oldNorm == newNorm )
			this.times(newNorm/oldNorm, excludePadding);
	}
	
	/**
	 * \brief Sets the norm to 1.0
	 * 
	 * @param excludePadding
	 */
	public void normalise(Boolean excludePadding)
	{
		this.normalise(1.0, excludePadding);
	}
	
	
	/*************************************************************************
	 * METHODS INVOLVING ANOTHER ARRAY
	 ************************************************************************/

	/**
	 * \brief Check if this array is the same as that given.
	 * 
	 * Will return false if any corresponding values in the arrays are different.
	 * 
	 * @param other	NDimensionalDouble object to be compared with this one.
	 * @return	Boolean stating whether other is the same (True) or different
	 * (False).
	 * @throws Exception if the dimensions are different.
	 */
	public Boolean equals(NDimensionalDouble other) throws Exception
	{
		/*
		 * Safety checking first.
		 */
		if ( ! areDimensionsSame(other) )
			throw new Exception();
		/*
		 * Check if any are different.
		 */
		for (int i = 0; i < this.totalLength; i++)
			if ( this.getValue(i) != other.getValue(i) )
				return false;
		/*
		 * If not, then the two arrays are the same.
		 */
		return true;
	}
	
	/**
	 * \brief 
	 * 
	 * @param other
	 * @throws Exception if the dimensions are different.
	 */
	public void add(NDimensionalDouble other, Boolean excludePadding)
															throws Exception
	{
		this.applyToAll((v1, v2) -> {return v1 + v2;}, other, excludePadding);
	}
	
	/**
	 * \brief 
	 * 
	 * @param other
	 * @throws Exception if the dimensions are different.
	 */
	public void subtract(NDimensionalDouble other, Boolean excludePadding)
															throws Exception
	{
		this.applyToAll((v1, v2) -> {return v1 - v2;}, other, excludePadding);
	}
	
	/**
	 * \brief 
	 * 
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception if the dimensions are different.
	 */
	public Double dotProduct(NDimensionalDouble other, Boolean excludePadding)
															throws Exception
	{
		return this.getProcessedBiSum((v1, v2) -> {return v1 * v2;},
													other, excludePadding);
	}
	
	/**
	 * \brief 
	 * 
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception if the dimensions are different.
	 */
	public Double distance(NDimensionalDouble other, Boolean excludePadding)
															throws Exception
	{
		Double sqSum = this.getProcessedBiSum(
								(v1, v2) -> {return ExtraMath.sq(v1 - v2);},
													other, excludePadding);
		return Math.sqrt(sqSum);
	}
	
	
	
	/*************************************************************************
	 * VECTOR METHODS
	 * These always assume there to be no padding.
	 ************************************************************************/
	
	/**
	 * \brief 
	 * 
	 * 
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception if the dimensions are different.
	 */
	public Double cosAngle(NDimensionalDouble other) throws Exception
	{
		Double dotProd = dotProduct(other, false);
		return ( dotProd == 0.0 ) ? 0.0 : dotProd/
										(this.norm(false) * other.norm(false));
	}
	
	/**
	 * \brief 
	 * 
	 * @param other
	 * @param excludePadding
	 * @return
	 * @throws Exception if the dimensions are different.
	 */
	public Double angle (NDimensionalDouble other) throws Exception
	{
		return Math.acos(this.cosAngle(other));
	}
	
	
}
