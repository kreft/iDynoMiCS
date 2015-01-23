package simulator.geometry.boundaryConditions;

import simulator.geometry.DiscreteVector;

public abstract class InternalBoundary extends AllBC
{
	public void applyBoundary(DiscreteVector coord)
	{
		/*
		 *  Nothing to do here: coordinates either side of an internal boundary
		 *  are valid.
		 */ 
	}
}
