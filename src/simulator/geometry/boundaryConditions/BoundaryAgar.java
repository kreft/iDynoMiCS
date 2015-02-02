package simulator.geometry.boundaryConditions;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import utils.XMLParser;

public class BoundaryAgar extends InternalBoundary
{
	public void init(Simulator aSim, Domain aDomain, XMLParser root)
	{
		readGeometry(root, aDomain);
		aDomain.addBoundary(this);
		_isSupport = true;
		
	}
	
	/**
	 * \brief If modeling an agar plate boundary, this method updates the boundary.
	 * 
	 * @param soluteGrid	Grid of all solutes
	 * @param reactionGrid	Grid of all reactions
	 * @param timeStep	Current internal timestep of the simulation
	 */
	public void updateAgar(SoluteGrid[] soluteGrid,
								SoluteGrid[] reactionGrid, Double timeStep)
	{
		
	}

	@Override
	public void setBoundary(LocatedGroup aGroup)
	{
		
		
	}

	@Override
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target)
	{
		hardBoundary(anAgent, target);
	};
}
