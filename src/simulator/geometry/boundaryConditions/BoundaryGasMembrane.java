/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  Project iDynoMicS
 * ___________________________________________________________________________
 * 
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * 
 */
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

import org.jdom.Element;
import java.util.*;
import utils.UnitConverter;
import utils.XMLParser;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;

/**
 * \brief BoundaryMembrane : defines a boundary impermeable to everything except to gas
 * 
 * BoundaryMembrane : defines a boundary impermeable to everything except to gas. A membrane boundary has a selective permeability,
 * meaning it behaves like a zero-flux boundary for agents and most of the solutes, but for selected solutes includes specification of 
 * the diffusivity in the membrane and the opposing-side solute concentration
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BoundaryGasMembrane extends AllBC
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long         serialVersionUID = 1L;

	/**
	 * The list of solutes to let diffuse through the membrane
	 */
	protected boolean[]               isPermeableTo;

	/**
	 * Level of permeability for each solute that can diffuse through the membrane
	 */
	protected double[]                permeability;

	/**
	 * The defined bulk in the simulation to which the liquid phase is connected
	 */
	protected Bulk                    _connectedBulk;

	/**
	 * A vector normal to the boundary and starting from the orthogonal projection
	 */
	protected static ContinuousVector vectorIn;

	/**
	 * \brief Declare a gas membrane boundary and set hasBulk to true to note this is the case
	 * 
	 * Declare a gas membrane boundary and set hasBulk to true to note this is the case
	 */
	public BoundaryGasMembrane() {
		hasBulk = true;
	}

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file, and builds the list of solutes to let diffuse through the membrane
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file, and builds the list of solutes to let diffuse through the membrane
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBCMarkUp	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBCMarkUp) {

		// this part is same as zero-flux boundary
		readGeometry(aBCMarkUp, aDomain);
		aDomain.addBoundary(this);
		_isSupport = true;

		// now need to set up the solute permeability
		
		String bulkName, soluteName;

		// Load description of the bulk connected to the membrane
		bulkName = aBCMarkUp.getParam("bulk");
		_connectedBulk = aSim.world.getBulk(bulkName);

		// Build the list of solutes to let diffuse through the membrane
		isPermeableTo = new boolean[aSim.soluteDic.size()];
		permeability = new double[aSim.soluteDic.size()];
		Arrays.fill(isPermeableTo, false);

		for (Element aChild : aBCMarkUp.getChildren("param")) {
			if (!aChild.getAttributeValue("name").equals("isPermeableTo")) continue;
			soluteName = aChild.getAttributeValue("detail");
			isPermeableTo[aSim.getSoluteIndex(soluteName)] = true;
			
			StringBuffer unit=new StringBuffer("");
			double paramValue = aBCMarkUp.getParamDbl("isPermeableTo", unit);
			paramValue *= UnitConverter.time(unit.toString());
			paramValue *= UnitConverter.length(unit.toString());
			paramValue *= UnitConverter.length(unit.toString());			
			permeability[aSim.getSoluteIndex(soluteName)] = paramValue;
		}
	}

	
	/**
	 * \brief Computes and applies gas diffusivity across the gas membrane boundary
	 * 
	 * Computes and applies gas diffusivity across the gas membrane boundary
	 * 
	 * @param relDif	Supplied RelDiff grid
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshDiffBoundary(SoluteGrid relDif, SoluteGrid aSoluteGrid) {
		double diffusivity;

		//Compute pseudo local diffusivity
		if (isPermeableTo[aSoluteGrid.soluteIndex]) {
			diffusivity = permeability[aSoluteGrid.soluteIndex]/aSoluteGrid.diffusivity;
		} else {
			diffusivity = 1;
		}
		
		// Apply or restore standard relative diffusivity
		_myShape.readyToFollowBoundary(relDif);
		while (_myShape.followBoundary(dcIn, dcOut, relDif)) {
			relDif.setValueAt(diffusivity, dcOut);
		}

	}

	/**
	 * \brief Solver for the gas membrane boundary condition. Initialises the course along the shape of the boundary. 
	 * 
	 * Solver for the gas membrane boundary condition. Initialises the course along the shape of the boundary
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) {

		// Initialise the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		if (isPermeableTo[aSoluteGrid.soluteIndex]) {
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
				aSoluteGrid.setValueAt(_connectedBulk.getValue(aSoluteGrid.soluteIndex), dcOut);
			}

		} else {
			// The membrane has the same behaviour than a zero-flux boundary
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
			}
		}
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
	public double getBulkValue(int soluteIndex) {
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
	public ContinuousVector lookAt(ContinuousVector cc) {
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
		aGroup.status = 0;
		// status 0 -> carrier
	}

	/**
     * Modify the movement vector : the new position is the orthognal projection
     * of the outside point on the boundary surface
     * 
     * @see LocatedAgent.move();
     */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) {
		// Define coordinates of the corrected position
		_myShape.orthoProj(target, target);

		// Build a vector normal to the boundary and starting from the
		// orthogonal projection
		vectorIn = new ContinuousVector(_myShape.getNormalInside(target));

		// The whole cell has to be inside, so make a translation equal to the
		// total radius
		vectorIn.times(anAgent.getRadius(true));

		// Compute the new position
		target.add(vectorIn);

		// Compute and update the movement vector leading to this new position
		anAgent.getMovement().sendDiff(anAgent.getLocation(), target);
	}

}
