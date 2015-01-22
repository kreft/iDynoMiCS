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

import simulator.Simulator;
import simulator.agent.LocatedGroup;
import utils.ExtraMath;

/**
 * \brief Quadratic detachment method. Function: kDet*L^2, where L is the
 * local biomass thickness.
 * 
 * The mark-up defines the erosion forces that act on the biofilm surface.
 * Detachment works by removing a layer of biomass based on the detachment
 * speed and the timestep, with the detachment speed calculated via one of the
 * given forms. This class captures the Quadratic detachment method.
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA).
 */
public class DS_Quadratic extends LevelSet
{
	/**
	 * \brief Calculate and return the local detachment speed using this
	 * detachment method.
	 *
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aGroup	Located group for which the local detachment speed is
	 * being determined.
	 * @return Double stating local detachment speed for this group.
	 */
	@Override
	protected Double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim)
	{
		if ( super.getLocalDetachmentSpeed(aGroup, aSim) == Double.MAX_VALUE )
			return Double.MAX_VALUE;
		return kDet*ExtraMath.sq(aGroup.cc.x);
	}
}