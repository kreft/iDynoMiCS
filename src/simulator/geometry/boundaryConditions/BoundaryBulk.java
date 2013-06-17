/**
 * \package simulator.geometry.boundaryConditions
 * \brief Package of boundary conditions that can be used to capture agent behaviour at the boundary of the computation domain
 * 
 * Package of boundary conditions that can be used to capture agent behaviour at the boundary of the computation domain. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry.boundaryConditions;

import utils.XMLParser;

import simulator.*;
import simulator.geometry.*;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;

/**
 * \brief Defines the bulk boundary: the concentration on the boundary is fixed by a dynamical bulk, the agents crossing this line die
 * 
 * Defines the bulk boundary: the concentration on the boundary is fixed by a dynamical bulk, the agents crossing this line die.
 * This boundary simulates the connection to a larger bulk liquid subject to a dilution process, as would be typical for a reactor. 
 * Behaviour for agents does not differ from the constant concentration boundary case (and so within the computational domain the 
 * boundary condition is identical), but solute dynamics in the bulk compartment require an additional computational step. In this step, 
 * ordinary differential equations describing the reactions occurring in the biofilm and the hydraulic processes affecting the bulk 
 * liquid are solved in order to determine the bulk concentration at the next time-step.
 *
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */
public class BoundaryBulk extends AllBC{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * The defined bulk in the simulation to which the liquid phase is connected
	 */
	protected Bulk            _connectedBulk;
	
	/**
	 * Value of solute in the bulk
	 */
	static double             bulkValue;

	/**
	 * \brief Declare a variable concentration boundary and set hasBulk to true to note this is the case
	 * 
	 * Declare a variable concentration boundary and set hasBulk to true to note this is the case
	 */
	public BoundaryBulk() {
		hasBulk = true;
	}

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file. In this case also links the connected bulk to this boundary
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file. In this case also links the connected bulk to this boundary
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBoundCondMarkUp	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBoundCondMarkUp) 
	{
		// Load the geometry of the boundary
		readGeometry(aBoundCondMarkUp, aDomain);		
		
		aDomain.addBoundary(this);
		
		// Load description of the connected bulk
		String bulkName = aBoundCondMarkUp.getParam("bulk");
		_connectedBulk = aSim.world.getBulk(bulkName);
	}
	
	/**
	 * \brief Solver for the variable concentration boundary condition. Initialises the course along the shape of the boundary. 
	 * 
	 * Solver for the variable concentration boundary condition. Initialises the course along the shape of the boundary
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) 
	{
		// Store the concentration in the bulk
		bulkValue = _connectedBulk.getValue(aSoluteGrid.soluteIndex);
		// bulkValue = 1;
		// Initialise the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
			aSoluteGrid.setValueAt(bulkValue, dcOut);
		}
	}

	/**
	 * \brief Updates the levels in the bulk. Allows reaction or flux-based bulk treatment
	 * 
	 *  Updates the levels in the bulk. Allows reaction or flux-based bulk treatment (BVM 151208)
	 * 
	 *  @param allSG	Array of all solute grids
	 *  @param allRG	Array of all reaction grids	
	 *  @param timeStep	The internal timestep currently being applied in this simulation
	 *  
	 */
	public void updateBulk(SoluteGrid[] allSG, SoluteGrid[] allRG, double timeStep) {
		_connectedBulk.updateBulk(allSG, allRG, timeStep);
	}

	/**
	 * \brief Return the bulk that is connected to this boundary
	 * 
	 * Return the bulk that is connected to this boundary
	 * 
	 * @return Bulk object that is connected to this boundary
	 */
	public Bulk getBulk() {
		return _connectedBulk;
	}

	/**
	 * \brief For a specified solute, returns the level of that solute in the bulk
	 * 
	 * For a specified solute, returns the level of that solute in the bulk
	 * 
	 * @param soluteIndex	Index of the solute in the simulation dictionary
	 * @return	Value of solute in the connected bulk
	 */
	public double getBulkValue(int soluteIndex) 
	{
		return _connectedBulk.getValue(soluteIndex);
	}

	/**
	 * \brief Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions need the input corrected, some don't and just return the input
	 * 
	 * Method used by another which gets the indexed grid position of a continuous vector. Some boundary conditions (e.g. BoundaryCyclic_ 
	 * need the input corrected due to the condition, some don't and just return the input. Maybe we'll change this at some point as to 
	 * just return the input looks a bit daft - but we'll leave it here for the moment
	 * 
	 * @param cc	ContinuousVector that gives the current location of an agent to check on the grid
	 */
	public ContinuousVector lookAt(ContinuousVector cc) 
	{
		return cc;
	}

	/**
     * \brief Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * @param aGroup	LocatedGroup object which has been detected to be outside the boundary
     */
	public void setBoundary(LocatedGroup aGroup) {
		aGroup.status = 3;
		// status 3 -> bulk
	}

	/**
	 * \brief Kills any agents that are crossing this boundary as they are leaving the simulated system
	 * 
	 * Kills any agents that are crossing this boundary as they are leaving the simulated system
	 * 
	 * @param anAgent	The LocatedAgent that has crossed the boundary
	 * @param target	Vector of where this agent was going to be placed
	 */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) 
	{
		//sonia 27.04.2010
		//recording reason of death (agent will be moved to agentToKill list when die() calls registerDeath()
		anAgent.death = "overBoard";
		
		anAgent.die(false);
		//LogFile.writeLog("agent killed by Bulk Boundary");
		
		// to label this agent as "shoving solved", set to zero its movement.
		anAgent.getMovement().reset();
		target.set(anAgent.getLocation());
		
	}


	/**
	 * \brief Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * @return String noting the side of the domain that this condition applies to (i.e. x0z, xNz, etc)
	 */
	public String toString() 
	{
		return new String("Bulk:"+this._mySide);
	}
}
