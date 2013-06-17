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

import java.util.List;

import org.jdom.Element;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;

import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import simulator.geometry.shape.IsShape;
import utils.XMLParser;

/**
 * \brief BoundaryCyclic : close the system along a dimension
 * 
 * BoundaryCyclic : close the system along a dimension. It is computationally unfeasible to simulate a micro-scale world on a 
 * macro-scale level, so oftentimes a small spatial sub-region is assumed to represent the system as a whole; in this case, 
 * periodic boundaries are used to remove artificial edge effects by assuming the simulated region adjoins other, similar, regions. 
 * As a consequence, boundaries in some chosen directions (generally for movements parallel to the substratum) are periodic, which 
 * means that the solute concentrations and solute gradients are constant across the boundary, and that agents travelling through one 
 * boundary will be translated to the other side of the domain
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BoundaryCyclic extends AllBC
{

	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long       serialVersionUID = 1L;

	/**
	 * Shape object containing a construction of the boundary opposite this one
	 */
	private IsShape                 _myOppShape;
	
	/**
	 * Vector that stores the intersection with the crossed boundary
	 */
	private static ContinuousVector vectorIn;
	
	/**
	 * Used to translate a set of points to their respective points on the opposite side of the boundary
	 */
	private static DiscreteVector   translator       = new DiscreteVector();

	/**
	 * \brief Declare a cyclic boundary and set isCyclic to true to note this is the case
	 * 
	 * Declare a cyclic boundary and set isCyclic to true to note this is the case
	 */
	public BoundaryCyclic() {
		isCyclic = true;
	}

	/**
	 * \brief Initialises the boundary from information contained in the simulation protocol file
	 * 
	 * Initialises the boundary from information contained in the simulation protocol file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aDomain	The domain which this boundary condition is associated with
	 * @param aBCParser	The XML tags that have declared this boundary in the protocol file
	 */
	public void init(Simulator aSim, Domain aDomain, XMLParser aBCParser) 
	{
		_mySide = aBCParser.getAttributeStr("name");
		
		// in 3D, all cyclic boundaries are active
		if(aDomain.is3D) activeForSolute=true;
		
		// in 2D, the x0y/xNy boundary is not active
		if(!aDomain.is3D & _mySide.contains("x0y")) activeForSolute=false;
		
		//activeForSolute = aBCParser.getParamBool("activeForSolute");
		readGeometry(aBCParser, aDomain);
		aDomain.addBoundary(this);
		aDomain.addBoundary(this.createOtherSide());		
	}

	/**
	 * \brief Read the geometry of this boundary from the protocol file and construct both this boundary and the opposite side (as this is cyclic)
	 * 
	 * Read the geometry of this boundary from the protocol file and construct both this boundary and the opposite side (as this is cyclic)
	 * 
	 * @see Domain constructor
	 * @param geometryRoot	XML tags in the protocol file that describe this boundary
	 * @param aDomain	The domain which this boundary is associated with
	 */
	public void readGeometry(XMLParser geometryRoot, Domain aDomain) {
		List<Element> shapeList = geometryRoot.getChildren("shape");
		String className;

		try {
			// Build first shape;
			className = "simulator.geometry.shape.";
			className += shapeList.get(0).getAttributeValue("class");
			_myShape = (IsShape) Class.forName(className).newInstance();
			_myShape.readShape(new XMLParser(shapeList.get(0)), aDomain);
			_mySide = geometryRoot.getAttributeStr("name");

			// Build opposite side shape

			className = "simulator.geometry.shape.";
			className += shapeList.get(1).getAttributeValue("class");
			_myOppShape = (IsShape) Class.forName(className).newInstance();
			_myOppShape.readShape(new XMLParser(shapeList.get(1)), aDomain);

		} catch (Exception e) {
		}
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
		ContinuousVector nCC = _myShape.intersection(cc, _myShape.getNormalInside(cc));
		ContinuousVector bCC = getSymmetric(nCC);
		bCC.subtract(nCC);
		nCC.sendSum(bCC, cc);
		return nCC;
	}

	/**
     * \brief Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * Change the status of a specified LocatedGroup to note that it has been identified as being outside this boundary
     * 
     * @param aGroup	LocatedGroup object which has been detected to be outside the boundary
     */
	public void setBoundary(LocatedGroup aGroup) {
		aGroup.status = -1;
		// status -1 -> outside
	}

	/**
	 * \brief Applies the boundary condition by modifying the movement vector. New position is orthogonal projection of the outside point on the boundary surface
	 * 
	 * Applies the boundary condition by modifying the movement vector. New position is orthogonal projection of the outside point on the boundary surface
	 * 
	 * @param anAgent	The Located Agent which is attempting to cross the boundary
	 * @param target	The target position that the agent is moving to
	 */
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) {
		// Determine the intersection with the crossed boundary
		vectorIn = _myShape.intersection(anAgent.getLocation(), anAgent.getMovement());

		// Determine the remaining movement when we touch the boundary
		target.sendDiff(target, vectorIn);

		// Apply the residual movement on the symmetric point
		vectorIn = getSymmetric(vectorIn);
		target.add(vectorIn);

		// Compute and update the movement vector leading to this new position
		anAgent.getMovement().sendDiff(target, anAgent.getLocation());
	}

	/**
	 * \brief Solver for the cyclic boundary condition. Initialises the course along the shape of the boundary, setting the values of solute near the boundary as required 
	 * 
	 * Solver for the cyclic boundary condition. Initialises the course along the shape of the boundary, , setting the values of solute near the boundary as required
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be refreshed by the solver
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) 
	{
		// 2D simulations: activeForSolute is false for x0y/xNy and true for x0z/xNz
		// 3D simulations: activeForSolute is always true for cyclic boundaries
		if (!activeForSolute) {
			// Initialise the course along the shape of the boundary
			_myShape.readyToFollowBoundary(aSoluteGrid);
			
			translator.set(_myOppShape.getNormalDC());
			translator.times(2);
			// Send a point belonging to the boundary and the closest point
			// outside the domain
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) 
			{
				
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
				DiscreteVector dcOutTemp = dcOut.clone();
				dcOut.add(translator);
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
			}
		} else {
			// Build translator between both boundaries
			int k = (int) Math.floor(_myOppShape.getDistance(_myShape)/aSoluteGrid.getResolution());
			translator.set(_myOppShape.getNormalDC());
			translator.times(k-1);

			// Initialise the course along the shape of the boundary
			_myShape.readyToFollowBoundary(aSoluteGrid);

			// Send a point belonging to the boundary and the closest point
			// outside the domain
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
				dcIn.add(translator);
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
			}
		}
	}

	/**
	 * \brief Applies boundary conditions for all other BC types - but not used for cyclic boundaries.
	 * 
	 * Applies boundary conditions for all other BC types - but not used for cyclic boundaries.
	 * 
	 * @param aSpGrid	SoluteGrid	
	 * @param dC	DiscreteVector
	 */
	public void applyBoundary(SoluteGrid aSpGrid, DiscreteVector dC) {
	}

	/**
	 * \brief Returns the value within a boundary. Not used for cyclic boundaries
	 * 
	 * Returns the value within a boundary. Not used for cyclic boundaries
	 * 
	 * @param aSpGrid	Solute grid which the boundary is being applied
	 * @param i	I Coordinate of grid space to query
	 * @param j	J Coordinate of grid space to query
	 * @param k	K Coordinate of grid space to query
	 * @return	Double balue of the solute level stored in that boundary cell
	 */
	public double getValueFromBoundary(SoluteGrid aSpGrid, int i, int j, int k) {
		System.out.println("CyclicBoundary:should not be used");
		return 0;
	}

	/**
	 * \brief Creates the opposite side of the cyclic boundary such that agents can 'roll around' to the other side
	 * 
	 * Creates the opposite side of the cyclic boundary such that agents can 'roll around' to the other side
	 * 
	 * @return	BoundaryCyclic object that captures the opposite side of this boundary
	 */
	public BoundaryCyclic createOtherSide() {
		BoundaryCyclic out = new BoundaryCyclic();
		out.activeForSolute = this.activeForSolute;
		out._myShape = this._myOppShape;
		out._myOppShape = this._myShape;
		out._mySide = this._mySide.replaceFirst("0", "N");
		return out;
	}

	/**
	 * \brief Determines if a point is outside the boundary
	 * 
	 * Determines if a point is outside the boundary
	 * 
	 * @param cc	ContinuousVector to check
	 * @return	Boolean value noting whether this coordinate is outside the boundary (true) or not (false)
	 */
	public boolean isOutside(ContinuousVector cc) {
		return _myShape.isOutside(cc);
	}

	/**
	 * \brief Determines if a point is on the boundary
	 * 
	 * Determines if a point is on the boundary
	 * 
	 * @param cC	ContinuousVector to check
	 * @param aSpatialGrid	Spatial grid to check
	 * @return	Boolean value noting whether this coordinate is on the boundary (true) or not (false)
	 */
	public boolean isOnBoundary(ContinuousVector cC, SpatialGrid aSpatialGrid) {
		return (_myShape.isOnBoundary(cC, aSpatialGrid.getResolution()));
	}

	/**
	 * \brief Returns the distance from a point to the boundary
	 * 
	 * Returns the distance from a point to the boundary
	 * 
	 * @param cc	The continuous vector of points to calculate how far the point is from the boundary
	 * @return	Double value stating the distance fromt the point to the boundary
	 */
	public double getDistance(ContinuousVector cc) {
		return _myShape.getDistance(cc);
	}

	/**
	 * \brief Return the intersection between the opposite shape and a provided point
	 * 
	 * Return the intersection between the opposite shape and a provided point
	 * 
	 * @param cc	A position on a boundary
	 * @return	ContinuousVector containing the intersection between the opposite shape and provided point
	 */
	public ContinuousVector getSymmetric(ContinuousVector cc) {
		// Determine on which shape you have to compute your future coordinates
		return _myOppShape.intersection(cc, _myShape.getNormalInside(cc));
	}

	
	/**
	 * \brief Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * @return String noting the side of the domain that this condition applies to (i.e. x0z, xNz, etc)
	 */
	public String toString() {
		return new String("Cyclic:"+this._mySide);
	}
}
