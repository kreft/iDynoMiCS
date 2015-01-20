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

import java.util.LinkedList;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.agent.zoo.EpithelialCell;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import utils.XMLParser;

/**
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public class BoundaryEpithelium extends InternalBoundary
{
	LinkedList<EpithelialCell> cells;
	
	public void init(Simulator aSim, Domain aDomain, XMLParser root)
	{
		readGeometry(root, aDomain);
		aDomain.addBoundary(this);
		
		
	}
	
	public void refreshBoundary(SoluteGrid aSoluteGrid)
	{
		
	}
	
	public ContinuousVector lookAt(ContinuousVector cc)
	{
		
		return null;
	}
	
	public void setBoundary(LocatedGroup aGroup)
	{
		
	}
	
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector newLoc)
	{
		
	}
	
}
