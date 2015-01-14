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
