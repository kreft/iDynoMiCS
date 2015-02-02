package simulator.geometry.boundaryConditions;

import simulator.SoluteGrid;
import simulator.geometry.DiscreteVector;

public abstract class InternalBoundary extends AllBC
{
	@Override
	public void refreshBoundary(SoluteGrid aSoluteGrid)
	{
		
	}
	
	@Override
	public void applyBoundary(DiscreteVector coord)
	{
		/*
		 *  Nothing to do here: coordinates either side of an internal boundary
		 *  are valid.
		 */ 
	}
}
