/**
 * \package simulator.geometry.boundaryConditions
 * \brief Package of boundary conditions that can be used to capture agent
 * behaviour at the boundary of the computation domain.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.geometry.boundaryConditions;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import utils.XMLParser;

public class ExternalBoundary extends AllBC
{
	
	@Override
	public void init(Simulator aSim, Domain aDomain, XMLParser root)
	{
		
	}
	
	@Override
	public void refreshBoundary(SoluteGrid aSoluteGrid)
	{
		
	}
	
	@Override
	public void setBoundary(LocatedGroup aGroup)
	{
		
	}
	
	@Override
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc)
	{
		
	}

}
