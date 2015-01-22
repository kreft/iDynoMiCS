/**
 * \package simulator.detachment
 * 
 * \brief Package of classes that capture detachment of agents from the
 * biofilm.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.detachment;

import simulator.agent.LocatedGroup;
import simulator.Simulator;
import simulator.SoluteGrid;

/**
 * \brief Loop through solutes to find the maximum gradient and use this as a
 * multiplier to calculate detachment rate.
 * 
 * The mark-up defines the erosion forces that act on the biofilm surface.
 * Detachment works by removing a layer of biomass based on the detachment
 * speed and the timestep, with the detachment speed calculated via one of the
 * given forms. This class captures the Proportional detachment method.
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham, UK.
 */
public class DS_SolGrad extends LevelSet 
{
	/**
	 *\brief Calculate and return the local detachment speed using this
	 *detachment method.
	 *
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aGroup	Located group for which the local detachment speed is
	 * being determined.
	 * @return Double stating local detachment speed for this group.
	 */
	protected Double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim) 
	{
		Double out = super.getLocalDetachmentSpeed(aGroup, aSim);
		if ( out == Double.MAX_VALUE )
			return out;
		/*
		 * If group below threshold, loop through solutes in order to find
		 * the maximum gradient at this point. Note that we don't have to
		 * worry about pressure (which is considered a solute in
		 * iDynoMiCS), because its bulk value should be 0.
		 * In 2D we count the y-direction twice.
		 */
		Double bulkValue;
		SoluteGrid solGrid;
		Double gradient;
		for (int solute = 0; solute < aSim.soluteDic.size(); solute++)
		{
			bulkValue = aSim.world.getMaxBulkValue(solute);
			if ( bulkValue > 0.0 )
			{
				solGrid = aSim.soluteList[solute];
				if ( aSim.is3D )
					gradient = solGrid.getGradient(aGroup.cc).norm();
				else
					gradient = solGrid.getGradient2D(aGroup.cc).norm();
				out = Math.max(out, kDet*gradient/bulkValue);
			}
		}
		return out;
	}
}