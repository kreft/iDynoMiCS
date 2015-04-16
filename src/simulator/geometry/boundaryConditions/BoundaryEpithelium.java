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
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import simulator.geometry.pointProcess.Site;
import utils.XMLParser;

/**
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public class BoundaryEpithelium extends InternalBoundary
{
	public LinkedList<Site> sites;
	
	public void init(Simulator aSim, Domain aDomain, XMLParser root)
	{
		readGeometry(root, aDomain);
		aDomain.addBoundary(this);
		_isSupport = true;
	}
	
	
	
	public void setBoundary(LocatedGroup aGroup)
	{
		
	}
	
	@Override
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target)
	{
		hardBoundary(anAgent, target);
	}
	
}
