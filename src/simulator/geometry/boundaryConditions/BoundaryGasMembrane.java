/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  Project iDynoMicS
 * ___________________________________________________________________________
 * BoundaryMembrane : defines a boundary impermeable to everything except to gas
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
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
public class BoundaryGasMembrane extends AllBC{

	/* ____________________________ FIELDS ________________________________ */
	// Serial version used for the serialisation of the class
	private static final long         serialVersionUID = 1L;

	// Defines permeability properties of the membrane
	protected boolean[]               isPermeableTo;
	protected double[]                permeability;

	// At which bulk the membrane is connected
	protected Bulk                    _connectedBulk;

	/* ______________ INTERNAL TEMPORARY VARIABLES ____________________ */
	protected static ContinuousVector vectorIn;

	/* ________________________ CONSTRUCTOR _______________________________ */
	public BoundaryGasMembrane() {
		hasBulk = true;
	}

	/**
     * Initialise a boundary condition simulating a selective membrane,
     * permeable to gas compounds
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

	/* ________________________ SOLVER _______________________________ */
	public void refreshDiffBoundary(SoluteGrid relDif, SoluteGrid aSoluteGrid) {
		double value;

		//Compute pseudo local diffusivity
		if (isPermeableTo[aSoluteGrid.soluteIndex]) {
			value = permeability[aSoluteGrid.soluteIndex]/aSoluteGrid.diffusivity;
		} else {
			value = 1;
		}
		
		// Apply or restore standard relative diffusivity
		_myShape.readyToFollowBoundary(relDif);
		while (_myShape.followBoundary(dcIn, dcOut, relDif)) {
			relDif.setValueAt(value, dcOut);
		}

	}

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

	public Bulk getBulk() {
		return _connectedBulk;
	}

	public double getBulkValue(int soluteIndex) {
		return _connectedBulk.getValue(soluteIndex);
	}


	/* _______________________ LOCATED AGENTS ______________________________ */
	public ContinuousVector lookAt(ContinuousVector cc) {
		return cc;
	}

	/**
     * Label a LocatedGroup which has been identified being outside this
     * boundary
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
