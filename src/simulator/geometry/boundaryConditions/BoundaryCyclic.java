/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________________________________
 * BoundaryCyclic : close the system along a dimension
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
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

public class BoundaryCyclic extends AllBC{

	// Serial version used for the serialisation of the class
	private static final long       serialVersionUID = 1L;

	private IsShape                 _myOppShape;
	private static ContinuousVector vectorIn;
	private static DiscreteVector   translator       = new DiscreteVector();

	/* ________________________ CONSTRUCTOR _______________________________ */
	public BoundaryCyclic() {
		isCyclic = true;
	}

	public void init(Simulator aSim, Domain aDomain, XMLParser aBCParser) {
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
	 * Achieve the construction of the object
	 * @see Domain constructor
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

	/* _______________________ LOCATED AGENTS ______________________________ */
	
	/**
	 * 
	 */
	public ContinuousVector lookAt(ContinuousVector cc) {
		ContinuousVector nCC = _myShape.intersection(cc, _myShape.getNormalInside(cc));
		ContinuousVector bCC = getSymmetric(nCC);
		bCC.subtract(nCC);
		nCC.sendSum(bCC, cc);
		return nCC;
	}

	/**
	 * Label a LocatedGroup which has been identified being outside this
	 * boundary
	 */
	public void setBoundary(LocatedGroup aGroup) {
		aGroup.status = -1;
		// status -1 -> outside
	}

	/**
	 * Modify the movement vector : the new position is the orthognal projection
	 * on the boundary surface
	 * @see LocatedAgent.move()
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

	/* ___________________________ SOLVER __________________________________ */
	/**
	 * Parse all the points included into the boundary and update their value
	 * according to the kind of boundary condition
	 * @see ComputationDomain.refreshBoundary()
	 */
	public void refreshBoundary(SoluteGrid aSoluteGrid) {
		// 2D simulations: activeForSolute is false for x0y/xNy and true for x0z/xNz
		// 3D simulations: activeForSolute is always true for cyclic boundaries
		if (!activeForSolute) {
			// Initialise the course along the shape of the boundary
			_myShape.readyToFollowBoundary(aSoluteGrid);
			
			translator.set(_myOppShape.getNormalDC());
			translator.times(2);
			// Send a point belonging to the boundary and the closest point
			// outside the domain
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid)) {
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
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
	 * Just for compatibility reason
	 */
	public void applyBoundary(SoluteGrid aSpGrid, DiscreteVector dC) {
	}

	public double getValueFromBoundary(SoluteGrid aSpGrid, int i, int j, int k) {
		System.out.println("CyclicBoundary:should not be used");
		return 0;
	}

	public void applyBoundary(SoluteGrid aSpGrid, DiscreteVector dC1, DiscreteVector dC2) {
		// dC2 is a point inside along the opposite border, dC1 is a point
		// outside along this border
		aSpGrid.setValueAt(aSpGrid.getValueAt(dC2), dC1);
	}

	/* _________________________ TOOLBOX ___________________________________ */

	public BoundaryCyclic createOtherSide() {
		BoundaryCyclic out = new BoundaryCyclic();
		out.activeForSolute = this.activeForSolute;
		out._myShape = this._myOppShape;
		out._myOppShape = this._myShape;
		out._mySide = this._mySide.replaceFirst("0", "N");
		return out;
	}

	/**
	 * If the cosinus of the angle defined by the the vector (ref this
	 * point),(normal vector) is positive, then the point is outside
	 */
	public boolean isOutside(ContinuousVector cc) {
		return _myShape.isOutside(cc);
	}

	public boolean isOnBoundary(ContinuousVector cC, SpatialGrid aSpatialGrid) {
		return (_myShape.isOnBoundary(cC, aSpatialGrid.getResolution()));
	}

	public ContinuousVector getIntersection(ContinuousVector position, ContinuousVector vector) {
		return _myShape.intersection(position, vector);
	}

	public ContinuousVector getOrthoProj(ContinuousVector cc) {
		return _myShape.getOrthoProj(cc);
	}

	public IsShape getClosestBoundary(ContinuousVector cc) {
		if (_myShape.getDistance(cc)<_myOppShape.getDistance(cc)) {
			return _myShape;
		} else {
			return _myOppShape;
		}
	}

	public double getDistance(ContinuousVector cc) {
		return _myShape.getDistance(cc);
	}

	/**
	 * Search the corresponding coordinates on the opposite boundary
	 * @param cc : a position on a boundary
	 * @param aShape : the shape on which the point is located
	 * @return
	 */
	public ContinuousVector getSymmetric(ContinuousVector cc) {
		// Determine on which shape you have to compute your future coordinates
		return _myOppShape.intersection(cc, _myShape.getNormalInside(cc));
	}

	/**
	 * @deprecated
	 * @param aSpGrid
	 */
	public void sendDCOnBothBorder(SpatialGrid aSpGrid) {

		// out1 = new ArrayList<DiscreteCoordinate>();
		// out2 = new ArrayList<DiscreteCoordinate>();
		int i, j, k;
		int nI = aSpGrid.getGridSizeI();
		int nJ = aSpGrid.getGridSizeJ();
		int nK = aSpGrid.getGridSizeK();
		double u[][][] = aSpGrid.getGrid();

		if (_mySide.equals("x0y")) {
			k = -1;
			for (i = -1; i<nI+1; i++) {
				for (j = -1; j<nJ+1; j++) {
					u[i+1][j+1][k+1] = u[i+1][j+1][nK];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(i+1,j+1,nK),i+1,j+1,k+1);
					// out1.add(new DiscreteCoordinate(i, j, k));
					// out2.add(new DiscreteCoordinate(i, j, nK-1));
				}
			}
		}

		if (_mySide.equals("xNy")) {
			k = 0;
			for (i = -1; i<nI+1; i++) {
				for (j = -1; j<nJ+1; j++) {
					u[i+1][j+1][nK+1] = u[i+1][j+1][k+1];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(i+1,j+1,k+1),i+1,j+1,nK);
					// out1.add(new DiscreteCoordinate(i, j, nK));
					// out2.add(new DiscreteCoordinate(i, j, k));
				}
			}
		}

		if (_mySide.equals("x0z")) {
			j = -1;
			for (i = -1; i<nI+1; i++) {
				for (k = -1; k<nK+1; k++) {
					u[i+1][j+1][k+1] = u[i+1][nJ][k+1];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(i+1,nJ,k+1),i+1,j+1,k+1);
					// out1.add(new DiscreteCoordinate(i, j, k));
					// out2.add(new DiscreteCoordinate(i, nJ-1, k));
				}
			}
		}

		if (_mySide.equals("xNz")) {
			j = 0;
			for (i = -1; i<nI+1; i++) {
				for (k = -1; k<nK+1; k++) {
					u[i+1][nJ+1][k+1] = u[i+1][j+1][k+1];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(i+1,j+1,k+1),i+1,nJ,k+1);
					// out2.add(new DiscreteCoordinate(i, j, k));
					// out1.add(new DiscreteCoordinate(i, nJ, k));
				}
			}
		}

		if (_mySide.equals("y0z")) {
			i = -1;
			for (j = -1; j<nJ+1; j++) {
				for (k = -1; k<nK+1; k++) {
					u[i+1][j+1][k+1] = u[nI][j+1][k+1];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(nI,j+1,k+1),i+1,j+1,k+1);
					// out1.add(new DiscreteCoordinate(i, j, k));
					// out2.add(new DiscreteCoordinate(nI-1, j, k));
				}
			}
		}

		if (_mySide.equals("yNz")) {
			i = 0;
			for (j = -1; j<nJ+1; j++) {
				for (k = -1; k<nK+1; k++) {
					u[nI+1][j+1][k+1] = u[i+1][j+1][k+1];
					// aSpGrid.setValueAt(aSpGrid.getValueAt(i+1,j+1,k+1),nI,nJ,k+1);
					// out2.add(new DiscreteCoordinate(i, j, k));
					// out1.add(new DiscreteCoordinate(nI, j, k));
				}
			}
		}
	}

	public String toString() {
		return new String("Cyclic:"+this._mySide);
	}
}
