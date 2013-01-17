/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ___________________________________________________________________________
 * ComputationDomain : default class to describe a Computation domain defined 
 * by a set of boundary conditions and a set of bulks
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */

package simulator.geometry;

import utils.ExtraMath;
import utils.XMLParser;

import java.util.*;

import simulator.*;
import simulator.geometry.boundaryConditions.*;
import simulator.agent.LocatedGroup;

public class Domain implements IsComputationDomain {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID    = 1L;

	public String             domainName;
	protected Simulator       currentSim;

	// Dimensions of the simulated environment in micrometer
	public double             length_X, length_Y, length_Z;
	public boolean            is3D;

	//sonia:chemostat 19.02.2010 changed visbility to public
	public int       _nI,_nJ, _nK;

	protected int             _i, _j, _k;

	public double             _resolution;
	protected SpatialGrid     _domainGrid;
	protected SoluteGrid      _boundaryLayer, _diffusivityGrid, _biomassGrid;

	// List of all boundaries defined on this computation domain
	private LinkedList<AllBC> _boundaryList       = new LinkedList<AllBC>();

	protected double          _dilationBand;
	protected double          specificArea;
	protected double          _biofilmDiffusivity = 1d;

	/* _______________________ CONSTRUCTOR ______________________________ */

	/**
	 * XML-based Constructor
	 * @param cdRoot
	 */
	public Domain(Simulator aSim, XMLParser cdRoot) {

		domainName = cdRoot.getAttribute("name");
		currentSim = aSim;

		is3D = (cdRoot.getChildAttrDbl("grid", "nDim")==3);
		_resolution = cdRoot.getParamLength("resolution");
			//sonia:
		// set the size of the grid to 1 instead of reading it from the protocol file
		if(Simulator.isChemostat){
			_nI =1;
			_nJ =1;
			_nK =1; 
				
		}else{
		
		_nI = (int) cdRoot.getChildAttrDbl("grid", "nI");
		_nJ = (int) cdRoot.getChildAttrDbl("grid", "nJ");
		if (is3D) _nK = (int) cdRoot.getChildAttrDbl("grid", "nK");
		else _nK = 1;
		
		}

		length_X = _nI*_resolution;
		length_Y = _nJ*_resolution;
		length_Z = _nK*_resolution;

		// Create and initialise the domain grid
		_domainGrid = new SoluteGrid(_nI, _nJ, _nK, _resolution, "domainGrid", this);

		// specific area is given in m2/m3
		specificArea = cdRoot.getParamDbl("specificArea");

		// Create the boundary layer and initialise it at "inside everywhere"
		_dilationBand = cdRoot.getParamLength("boundaryLayer");
		_boundaryLayer = createGrid( "boundaryLayer", 1);
		
		// Create the biomass MASS grid and initialise it empty
		_biomassGrid = createGrid( "totalBiomass", 0);

		// Create the relative diffusivity grid and initialise it at "liquid
		// everywhere"
		_biofilmDiffusivity = cdRoot.getParamDbl("biofilmDiffusivity");
		_diffusivityGrid = createGrid( "diffusivityGrid", 1);

		// Create all boundary conditions
		for (XMLParser aBCMarkUp : cdRoot.buildSetParser("boundaryCondition"))
			AllBC.staticBuilder(aBCMarkUp, aSim, this);
	//sonia:chemostat
		if(Simulator.isChemostat){
			//do nothing...
		}else{
		// Build the domain grid : 0 outside, 1 inside, -1 carrier
		applyAllBoundary();
		}
	}

	public SoluteGrid createGrid(String gridName,double defaultValue){
		SoluteGrid aGrid = new SoluteGrid(_nI, _nJ, _nK, _resolution,gridName, this);
		aGrid.setAllValueAt(defaultValue);
		return aGrid;
	}
	
	/**
	 * To use when you have finished to register all boundary conditions ; this
	 * function will populate a 3-D matrix. Apply all boundaries one after one ;
	 * a point is outside of the computational domain if it is declared outside
	 * by at least one of the Boundary Conditions
	 */
	public void applyAllBoundary() {
		DiscreteVector dC = new DiscreteVector();
		ContinuousVector cC;

		// Reset all the computational domain to "inside";
		_domainGrid.setAllValueAt(1);

		for (int i = 0; i<_domainGrid.getGridTotalSize(1); i++) {
			for (int j = 0; j<_domainGrid.getGridTotalSize(2); j++) {
				repet: for (int k = 0; k<_domainGrid.getGridTotalSize(3); k++) {
					dC.set(i-1, j-1, k-1);
					cC = _domainGrid.getContinuousCoordinates(dC);

					for (AllBC aBC : _boundaryList) {
						// skip if this gridCell has already been updated
						if (_domainGrid.getValueAt(i, j, k)==-1) continue repet;

						// Test if this grid cell is seen outside
						if (aBC.isOutside(dC, _domainGrid)) {
							_domainGrid.setValueAt(-1, i, j, k);
							continue repet;
						}

						// label carrier part of the domain
						if (aBC.isSupport() && aBC.getDistance(cC)<_resolution)
							_domainGrid.setValueAt(0, i, j, k);
					}
				}
			}
		}
	}

	/* _______________________ LOCATED AGENTS ______________________________ */
	/**
     * Test if a given location is outside a boundary ; used to detect the
     * crossed boundary when moving an agent
	 * @param newLoc
	 * @see LocatetAgent.move()
     * @param aBoundary : reference used to store the crossed boundary
     * @return true if a boundary has been crossed
     */
	public AllBC testCrossedBoundary(ContinuousVector newLoc) {
		// Test on the domain grid if the new location is inside the domain
		if (_domainGrid.isValid(newLoc)&&_domainGrid.getValueAt(newLoc)>=0) return null;

		// Find the first of the boundaries which has been crossed
		for (int iB = 0; iB<_boundaryList.size(); iB++) {
			//System.out.println("\tiB: "+iB);
			//if (!_boundaryList.get(iB).isActive()) continue;
			if (_boundaryList.get(iB).isOutside(newLoc)) { return _boundaryList.get(iB); }
		}

		// If you are here, it means that no boundary is being crossed.
		return null;
	}

	/* __________________________ TOOLBOX _________________________________ */

	public boolean isInside(ContinuousVector cc) {
		for (AllBC aBC : _boundaryList) {
			if (aBC.isOutside(cc)) return false;
		}
		return true;
	}

	public void addBoundary(AllBC aBC) {
		_boundaryList.add(aBC);
	}

	public LinkedList<AllBC> getAllBoundaries() {
		return _boundaryList;
	}

	/**
	 * Refresh relative diffusivity and boundary layer grids
	 */
	public void refreshBioFilmGrids() {
		//sonia:chemostat
		if(Simulator.isChemostat){
		// Build a grid with the concentration of agents
		// skip the the refreshment of the position of the agents relative to the boundary layers	
		_biomassGrid.setAllValueAt(0);
		currentSim.agentGrid.fitAgentMassOnGrid(_biomassGrid);
		
		}else{
			
		
		// Build a grid with the concentration of agents
		_biomassGrid.setAllValueAt(0);
		currentSim.agentGrid.fitAgentMassOnGrid(_biomassGrid);

		// reset the grid
		_boundaryLayer.setAllValueAt(0d);

		for (int i = 1; i<_nI+1; i++) {
			for (int j = 1; j<_nJ+1; j++) {
				for (int k = 1; k<_nK+1; k++) {
					if (_biomassGrid.grid[i][j][k]>0) {
						// if this is biomass,
						_boundaryLayer.grid[i][j][k] = 1.0d;
						_diffusivityGrid.grid[i][j][k] = _biofilmDiffusivity;
					} else {
						// if liquid, check dilation sphere for biomass
						// (checkDilationRadius will set the value to 1 if it is
						//  within the boundary layer)
						_boundaryLayer.grid[i][j][k] = checkDilationRadius(i, j, k);
						if (_domainGrid.grid[i][j][k]==-1) {
							_diffusivityGrid.grid[i][j][k] = Double.MIN_VALUE;
						} else {
							_diffusivityGrid.grid[i][j][k] = 1.0d;
						}
					}
				}
			}
		}

		_boundaryLayer.refreshBoundary();
		_diffusivityGrid.refreshBoundary();
		_biomassGrid.refreshBoundary();
		}
	}

	/**
	 * @return a list of DC with the limit of the boundary layer
	 */
	public LinkedList<DiscreteVector> getBorder() {
		double v;
		LinkedList<DiscreteVector> border = new LinkedList<DiscreteVector>();
		for (_i = 1; _i<_nI+1; _i++) {
			for (_j = 1; _j<_nJ+1; _j++) {
				for (_k = 1; _k<_nK+1; _k++) {
					v = _boundaryLayer.grid[_i][_j][_k];
					if ((v==1) && bdryHasFreeNbh()) {
						// add the location if it has biomass or is in the boundary layer (v==1) and
						// if the neighboring points are free (not biomass or bdry layer)
						border.addLast(new DiscreteVector(_i, _j, _k));
					}
				}
			}
		}
		return border;
	}

	/**
	 * @return a list of Doubles with the heights of the biofilm/liquid interface
	 * 			(this routine is very basic in that it just captures the overall
	 * 			 height (x-position) of each point from the nearest carrier,
	 * 			 but for most cases it should be sufficient)
	 */
	public double[] getInterface() {
		currentSim.agentGrid.getLevelSet().refreshBorder(false, currentSim);
		LinkedList<LocatedGroup> border =
				currentSim.agentGrid.getLevelSet().getBorder();
		
		
		// catch if there is no biomass for some reason; in that case return zero height
		if (border.size() == 0)	return new double[]{0.};

		// now copy to regular array, but watch for infinite distances
		double [] borderarray = new double[border.size()];
		for (int i = 0; i < border.size(); i++) {
			borderarray[i] = border.get(i).distanceFromCarrier;
			if (borderarray[i] == Double.MAX_VALUE) borderarray[i] = 0.;
		}

		return borderarray;
	}

	// bvm note 16.12.08: renamed this from hasFreeNbh() to differentiate from the 
	// biofilm version added below
	private boolean bdryHasFreeNbh() {
		if (is3D()) {
			// 3D grid
			for (int i = -1; i<2; i++) {
				for (int j = -1; j<2; j++) {
					for (int k = -1; k<2; k++) {
						if (_boundaryLayer.grid[_i+i][_j+j][_k+k]==0) return true;
					}
				}
			}
			return false;
		} else {
			// 2D grid
			for (int i = -1; i<2; i++) {
				for (int j = -1; j<2; j++) {
					if (_boundaryLayer.grid[_i+i][_j+j][1]==0) return true;
				}
			}
			return false;
		}
	}

	// bvm note 16.12.08: copied this routine following above for the boundary
	// returns true only if there is a free neighbor in the biomass grid
	private boolean bflmHasFreeNbh() {
		if (is3D()) {
			// 3D grid
			for (int i = -1; i<2; i++) {
				for (int j = -1; j<2; j++) {
					for (int k = -1; k<2; k++) {
						if (_biomassGrid.grid[_i+i][_j+j][_k+k] <= 0) return true;
					}
				}
			}
			return false;
		} else {
			// 2D grid
			for (int i = -1; i<2; i++) {
				for (int j = -1; j<2; j++) {
					if (_biomassGrid.grid[_i+i][_j+j][1] <= 0) return true;
				}
			}
			return false;
		}
	}
	
	protected int checkDilationRadius(int n, int m, int l) {
		// for no boundary layer, liquid means it's outside the boundary
		// (and this routine only checks the status of non-biomass elements)
		if (_dilationBand==0.) return 0;
		
		int nInterval, mInterval, lInterval;
		double deltaN, deltaM;
		double dilationRadiusM, dilationRadiusL;

		nInterval = (int) Math.floor(_dilationBand/_resolution);

		for (int i = (-nInterval); i<=nInterval; i++) {
			// only proceed if neighbour is within computational
			// volume top and bottom boundaries
			if ((n+i>=0)&(n+i<_nI)) {
				deltaN = i*_resolution;

				// This calculates the range in the j direction based on a right triangle
				// with hypotenuse equal to the sphere's radius, so that the total area
				// checked is a sphere
				dilationRadiusM = Math.sqrt(ExtraMath.sq(_dilationBand)-ExtraMath.sq(deltaN));
				mInterval = (int) Math.floor(dilationRadiusM/_resolution);

				for (int j = (-mInterval); j<=mInterval; j++) {
					if (_nK==1) {
						// 2D case
						if (_biomassGrid.grid[n+i][cyclicIndex(m+j,_nJ+2)][1]>0) { return 1; }
						if (_domainGrid.grid[n+i][cyclicIndex(m+j,_nJ+2)][1]==0) { return 1; }
					} else {
						// 3D case
						deltaM = j*_resolution;

						// This calculates the range in the k direction based on
						// a right triangle with hypotenuse equal to the sphere's
						// radius, so that the total area checked is a sphere
						dilationRadiusL = Math.sqrt(ExtraMath.sq(_dilationBand)-
									ExtraMath.sq(deltaN)-ExtraMath.sq(deltaM));
						lInterval = (int) Math.floor(dilationRadiusL/_resolution);

						for (int k = (-lInterval); k<=lInterval; k++) {
							if ((i!=0)|(j!=0)|(k!=0)) {
								if (_biomassGrid.grid[n+i]
								                     [cyclicIndex(m+j,_nJ+2)]
								                     [cyclicIndex(l+k,_nK+2)]>0)
									return 1;
								if (_domainGrid.grid[n+i]
								                    [cyclicIndex(m+j,_nJ+2)]
								                    [cyclicIndex(l+k,_nK+2)]==0)
									return 1;
							}
						}
					}
				}
			}
		}
		return 0;
	}

	protected final int cyclicIndex(int val, int limit) {
		return (val<0 ? limit+val : (val>=limit ? val-limit : val));
	}

	/* _____________________ ACCCESSOR ___________________________ */
	public double getLongestSize() {
		return ExtraMath.max(ExtraMath.max(length_X, length_Y), length_Z);
	}

	public double getResolution() {
		return _resolution;
	}

	/**
	 * @return the domain grid
	 */
	public SpatialGrid getGrid() {
		return _domainGrid;
	}

	public String getName() {
		return domainName;
	}

	public boolean is3D() {
		return _domainGrid.is3D();
	}

	public SoluteGrid getDiffusivity() {
		return _diffusivityGrid;
	}

	public SoluteGrid getBoundaryLayer() {
		return _boundaryLayer;
	}

	public SoluteGrid getBiomass() {
		return _biomassGrid;
	}

	/* ________________ UNUSED FUNCTIONS _________________________ */
	
	/**
	 * @deprecated
	 */
	public AllBC getClosestBoundary(ContinuousVector cc) {
		AllBC closeBC = _boundaryList.get(0);
		double d = closeBC.getDistance(cc);
		for (AllBC aBC : _boundaryList) {
			if (aBC.getDistance(cc)<d) {
				d = aBC.getDistance(cc);
				closeBC = aBC;
			}
		}
		return closeBC;
	}

	/**
	 * Return the closest boundary condition with a bulk Usefull to set
	 * concentration outside the diffusive domain
	 * @param cc
	 * @return
	 * 
	 * @deprecated
	 */
	public AllBC getClosestBoundaryWithBulk(ContinuousVector cc) {
		AllBC closeBC = _boundaryList.get(0);
		double d = 0;

		for (AllBC aBC : _boundaryList) {
			if (aBC.hasBulk()) {
				closeBC = aBC;
				d = closeBC.getDistance(cc);
				break;
			}
			return null;
		}

		for (AllBC aBC : _boundaryList) {
			if (aBC.hasBulk()&&aBC.getDistance(cc)<d) {
				d = aBC.getDistance(cc);
				closeBC = aBC;
			}
		}
		return closeBC;
	}
}
