/**
 * \package simulator.geometry
 * \brief Package of boundary utilities that aid the creation of the environment being simulated
 * 
 * Package of boundary utilities that aid the creation of the environment being simulated. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry;

import utils.ExtraMath;
import utils.XMLParser;

import java.util.*;

import simulator.*;
import simulator.geometry.boundaryConditions.*;
import simulator.agent.LocatedGroup;

/**
 * \brief Define the computation domain: an evenly spaced rectilinear grid described by its dimensionality (2D or 3D)
 * 
 * The computation domain is an evenly spaced rectilinear grid described by its dimensionality (2D or 3D), its size, geometry 
 * and the behaviour at its boundaries. See Figure 1 of the Lardon et al paper (2011) for a good description on how this is 
 * divided into several regions - the support, the bulk, the biofilm matrix, and the diffusion boundary layer. 
 * 
 * @since June 2006
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class Domain implements IsComputationDomain 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID    = 1L;

	/**
	 * Name of this computation domain, as supplied in the specified protocol file
	 */
	public String	domainName;
	
	/**
	 * The simulation object being used to simulate the conditions specified in the protocol file
	 */
	protected Simulator       currentSim;
	
	/**
	 * Domain X dimension in micrometers
	 */
	public double	length_X;
	
	/**
	 * Domain Y dimension in micrometers
	 */
	public double	length_Y;
	
	/**
	 * Domain Z dimension in micrometers
	 */
	public double	length_Z;
	
	/**
	 * Whether this computation domain is two or three dimensional
	 */
	public boolean	is3D;

	/**
	 * Number of grid elements in the x direction
	 */
	public int	_nI;
	
	/**
	 * Number of grid elements in the y direction
	 */
	public int	_nJ;
	
	/**
	 * Number of grid elements in the z direction
	 */
	public int 	_nK;

	protected int	_i;
	
	protected int	_j;
	
	protected int	_k;

	/**
	 * Width of each side of the grid element (in micrometres) 
	 */
	public double             _resolution;
	
	/**
	 * The solute grid that is a component of this computation domain
	 */
	protected SpatialGrid     _domainGrid;
	
	/**
	 * Boundary layer between bulk and biofilm
	 */
	protected SoluteGrid      _boundaryLayer;
	
	/**
	 * 	Diffusivity of solutes in each area of this domain
	 */
	protected SoluteGrid	_diffusivityGrid;
	
	/**
	 * Array to hold the X position that is the top of the boundary layer in each of the Y and Z positions. Used for self-attach scenarios
	 */
	public double[][]		_topOfBoundaryLayer;
	
	/**
	 * Grid to hold total biomass in a particular area
	 */
	public SoluteGrid 	_biomassGrid;

	/**
	 * List of all boundaries defined on this computation domain
	 */
	private LinkedList<AllBC> _boundaryList       = new LinkedList<AllBC>();
	
	/**
	 * Band between the boundary and bulk, capturing change in diffusivity and solute levels
	 */
	protected double          _dilationBand;
	
	/**
	 * The ratio between the carrier surface (the substratum on which the biofilm grows) and the bulk compartment volume. The physical
	 * volume of the system that appears in the simulation definition. In m2/m3.
	 */
	protected double          specificArea;
	
	/**
	 * Factor used to decrease solute diffusivity inside the biofilm. Multiplicative factor applied to the diffusivity in water
	 */
	protected double          _biofilmDiffusivity = 1d;

	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/

	/**
	 * \brief Creates a computation domain compartment object with attributes specified in the protocol file
	 * 
	 * The computation domain is an evenly spaced rectilinear grid described by its dimensionality (2D or 3D), its size, geometry 
	 * and the behaviour at its boundaries. See Figure 1 of the Lardon et al paper (2011) for a good description on how this is 
	 * divided into several regions - the support, the bulk, the biofilm matrix, and the diffusion boundary layer. The definition within 
	 * the computationDomain markup of the protocol file notes how these regions are set up. This constructor sets up each computation 
	 * domain that is specified in the protocol file.
	 * 
	 * @param aSim	The simulation object upon which the scenario specified in the protocol file is being run
	 * @param cdRoot	The XML tag objects that sub-nodes of the 'Bulk' tag in the protocol file
	 */
	public Domain(Simulator aSim, XMLParser cdRoot) 
	{
		domainName = cdRoot.getAttribute("name");
		currentSim = aSim;

		// Now determine if this computation domain is 2D or 3D
		is3D = (cdRoot.getChildAttrDbl("grid", "nDim")==3);
		_resolution = cdRoot.getParamLength("resolution");
		
		// Determine the size of this domain. If this is a chemostat, this is a set size and does not have to be read in.
		// If not a chemostat, take from XML file. (Added by Sonia Martins)
		if(Simulator.isChemostat)
		{
			_nI =1;
			_nJ =1;
			_nK =1; 
				
		}
		else
		{
			_nI = cdRoot.getChildAttrInt("grid", "nI");
			_nJ = cdRoot.getChildAttrInt("grid", "nJ");
			_nK = (is3D) ? cdRoot.getChildAttrInt("grid", "nK") : 1;
		}

		// Now calculate the length of the grid in micrometres
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
		
		//printBoundaryLayer();
		
		// Create the biomass MASS grid and initialise it empty
		_biomassGrid = createGrid( "totalBiomass", 0);

		// Create the relative diffusivity grid and initialise it at "liquid everywhere"
		_biofilmDiffusivity = cdRoot.getParamDbl("biofilmDiffusivity");
		_diffusivityGrid = createGrid( "diffusivityGrid", 1);

		// Create the grid that is going to hold the top of the boundary layer
		// Note that the boundary layer is initialised with padding, then ignores it (and never calculates the padding)
		// For consistency (for the moment), the top of the boundary layer array will work in the same way - thus we will 
		// be ignoring space at either end (_nJ and _nK=0 and _nJ+2 and _nK+2)  
		_topOfBoundaryLayer = new double[_nJ+2][_nK+2];
		
		// Now comes the definition of the behavior at the boundaries. In general, there are 6
		// boundaries that must be addressed: y0z, yNz, x0z, xNz, x0y, xNy. These represent
		// the edges of the domain along the non-named direction (i.e. y0z is the face at x=0,
		// and yNz is the face at x=N). (For 2D simulations the x0y and xNy directions are
		// included, but are made periodic.)
		// Each <boundaryCondition> also includes a <shape> mark-up to define the shape of the boundary.
		// The below call combines all boundary conditions in the XML file, then processes each
		for (XMLParser aBCMarkUp : cdRoot.buildSetParser("boundaryCondition"))
			AllBC.staticBuilder(aBCMarkUp, aSim, this);
		
		// Note the above has added all the boundaries to the array _boundaryList
	
		// Now apply these boundaries
		
		if(!Simulator.isChemostat)
		{
			// Build the domain grid : 0 outside, 1 inside, -1 carrier
			applyAllBoundary();
		}
		
		// KA May 2013
		// Now we're going to initialise all these grids
		// Note this wasn't previously done, but with self attachment we need to know where the boundary layer is before any agents
		// are added. The function below was part of the refreshBioFilmGrids method - this now exists independently and is called whenever 
		// these grids need to be refreshed
		calculateComputationDomainGrids();
		
		// Now we can initialise the top of the boundary layer
		this.calculateTopOfBoundaryLayer();
		
	}

	/**
	 * \brief Creates a solute or species grid and initialises the values within that grid. Used to create boundary and biomass grids
	 * 
	 * Creates a solute or species grid and initialises the values within that grid. Used to create boundary and biomass grids
	 * 
	 * @param gridName	The name of the grid being created (e.g. boundaryLayer / totalBioMass)
	 * @param defaultValue	The default value to assign to all grid spaces
	 * @return	Initialised solute grid of the size required by the simulation, initialised to the given default value
	 */
	public SoluteGrid createGrid(String gridName,double defaultValue)
	{
		SoluteGrid aGrid = new SoluteGrid(_nI, _nJ, _nK, _resolution,gridName, this);
		aGrid.setAllValueAt(defaultValue);
		return aGrid;
	}
	
	/**
	 * \brief Applies all specified domain boundaries to this computation domain one by one
	 * 
	 * This method should be used when you have finished to register all boundary conditions; this function will populate a 3-D matrix. 
	 * Apply all boundaries one after one ; a point is outside of the computational domain if it is declared outside by at least one 
	 * of the Boundary Conditions
	 * 
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
     * \brief Test if a given location is outside a boundary. Used to detect the crossed boundary when moving an agent
     * 
     * Test if a given location is outside a boundary. Used to detect the crossed boundary when moving an agent
     * 
     * @param newLoc	The location to test
     * @return	Boundary that the point has crossed (if applicable - null if no boundary crossed)
     */
	public AllBC testCrossedBoundary(ContinuousVector newLoc) {
		// Test on the domain grid if the new location is inside the domain
		if (_domainGrid.isValid(newLoc)&&_domainGrid.getValueAt(newLoc)>=0) return null;

		// Find the first of the boundaries which has been crossed
		for (int iB = 0; iB<_boundaryList.size(); iB++) 
		{
			if (_boundaryList.get(iB).isOutside(newLoc)) 
			{ 
				return _boundaryList.get(iB); 
			}
		}

		// If you are here, it means that no boundary is being crossed.
		return null;
	}
	
	/**
     * \brief Test if a given location is outside a boundary other than the bottom. Used in self-attach scenarios
     * 
     * Test if a given location is outside a boundary other than the bottom. For self-attachment, the simulation detects a swimming agent 
     * may have crossed the substratum boundary and will then assume that agent attaches. However we need to check that move has not crossed
     * any of the other boundaries, else that move is invalid. To do this, all boundaries are checked. If using the method above, y0z could 
     * still be returned and thus we end up in a loop. Thus this has had to be adapted so this cannot be returned
     * 
	 * @param newLoc	The new location of this swimming agent
     * @return The boundary which this move crosses, if applicable. Null if no such boundary
     */
	public AllBC testCrossedBoundarySelfAttach(ContinuousVector newLoc) 
	{
		// Test on the domain grid if the new location is inside the domain
		if (_domainGrid.isValid(newLoc)&&_domainGrid.getValueAt(newLoc)>=0) return null;

		// Find the first of the boundaries which has been crossed
		for (int iB = 0; iB<_boundaryList.size(); iB++) 
		{
			// Added a check to not return y0z - we know this has been crossed as the cell has met the substratum. We are only interested
			// here in checking the other 7 boundaries
			if (_boundaryList.get(iB).isOutside(newLoc)) 
			{ 
				if(_boundaryList.get(iB).getSideName() != "y0z");
				{
					return _boundaryList.get(iB);
				}
				 
			}
		}

		// If you are here, it means that no boundary is being crossed.
		return null;
	}

	/**
	 * \brief Add a boundary condition to the list of boundaries on this domain
	 * 
	 * Add a boundary condition to the list of boundaries on this domain
	 * 
	 * @param aBC	Boundary condition to add to the list of boundaries
	 */
	public void addBoundary(AllBC aBC) 
	{
		_boundaryList.add(aBC);
	}

	/**
	 * \brief Return all the boundary conditions applicable to this domain
	 * 
	 * Return all the boundary conditions applicable to this domain
	 * 
	 * @return LinkedList of boundary conditions for this domain
	 */
	public LinkedList<AllBC> getAllBoundaries() 
	{
		return _boundaryList;
	}
	
	/**
	 * \brief Creates an array list containing the 'i' coordinate of the computation domain that is the top of the boundary layer
	 * 
	 * Creates an array list containing the 'i' coordinate of the computation domain that is the top of the boundary layer. Note that 
	 * this is in the resolution specified by the 'resolution' parameter in the computation domain section of the protocol file, and this 
	 * may need to be adjusted if being used in calculations where the resolution differs
	 */
	public void calculateTopOfBoundaryLayer()
	{
		for(int k=1;k<_boundaryLayer.getGridSizeK()+1;k++)
		{
			for(int j=1;j<_boundaryLayer.getGridSizeJ()+1;j++)
			{
				int i=1;
				while(i<_boundaryLayer.getGridSizeI()+1 && _boundaryLayer.getValueAt(i, j, k)>0)
				{
					i++;
				}
				
				// assume now we've reached the point where 'i' has become 0, and thus we are out of the boundary layer
				// Subtract 1 such that the top of the layer is noted, not the outside of the layer
				_topOfBoundaryLayer[j][k]=i-1;
			}
		}
	}

	/**
	 * \brief Refresh relative diffusivity and boundary layer grids to ensure biomass updated this step is included
	 * 
	 * Method to refresh relative diffusivity and boundary layer grids to ensure biomass updated this step is included. 
	 * Used in the creation of output files
	 */
	public void refreshBioFilmGrids() 
	{
		//sonia:chemostat
		if(Simulator.isChemostat)
		{
			// Build a grid with the concentration of agents
			// skip the the refreshment of the position of the agents relative to the boundary layers	
			_biomassGrid.setAllValueAt(0);
			currentSim.agentGrid.fitAgentMassOnGrid(_biomassGrid);
			
			
		
		}
		else
		{
			// Build a grid with the concentration of agents
			_biomassGrid.setAllValueAt(0);
			
			currentSim.agentGrid.fitAgentMassOnGrid(_biomassGrid);
	
			// reset the grid
			_boundaryLayer.setAllValueAt(0d);
	
			// calculate the values in each of the grids
			calculateComputationDomainGrids();
			
			_boundaryLayer.refreshBoundary();
			
			// Now calculate the positions that are at the top of the boundary layer
			calculateTopOfBoundaryLayer();
				
			_diffusivityGrid.refreshBoundary();
			_biomassGrid.refreshBoundary();
			
		}
	}
	
	/**
	 * \brief Calculates the diffusivity and boundary layer grid levels
	 * 
	 * Calculates the diffusivity and boundary layer grid levels. In previous versions of iDynoMiCS this method could 
	 * be found within refreshBioFilmGrids. This has been moved here as, with the addition of self attachment, this method needs to be 
	 * called before agent initialisation. KA May 2013
	 */
	public void calculateComputationDomainGrids()
	{
		for (int i = 1; i<_nI+1; i++) 
		{
			for (int j = 1; j<_nJ+1; j++) 
			{
				for (int k = 1; k<_nK+1; k++) 
				{
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
	}

	
	/**
	 * \brief Returns a list of discrete vectors that specify the limit of the boundary layer
	 * 
	 * Returns a list of discrete vectors that specify the limit of the boundary layer
	 * 
	 * @return LinkedList of DiscreteVectors with the limit of the boundary layer
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
	 * \brief Creates a list of doubles with the heights of the biofilm/liquid interface. Used for writing simulation statistics
	 * 
	 * Creates a list of doubles with the heights of the biofilm/liquid interface. Used for writing simulation statistics. This routine 
	 * is very basic in that it just captures the overall height (x-position) of each point from the nearest carrier, but for most 
	 * cases it should be sufficient)
	 * 
	 * @return	Array of double values of the heights of the biofilm/liquid interface
	 */
	public Double[] getInterface()
	{
		currentSim.agentGrid.getLevelSet().refreshBorder(false, currentSim);
		LinkedList<LocatedGroup> border =
				currentSim.agentGrid.getLevelSet().getBorder();
		
		// Catch if there is no biomass for some reason; in that case return zero height
		if (border.size() == 0)	return ExtraMath.newDoubleArray(1);
		
		// Now copy to regular array, but watch for infinite distances
		Double [] borderarray = ExtraMath.newDoubleArray(border.size());
		for (int i = 0; i < border.size(); i++)
			if (border.get(i).distanceFromCarrier < Double.MAX_VALUE)
				borderarray[i] = border.get(i).distanceFromCarrier;
		
		return borderarray;
	}

	/**
	 * \brief Determines whether points in the boundary layer have free neighbours
	 * 
	 * Determines whether points in the boundary layer have free neighbours
	 * 
	 * @author BVM 161208
	 * @return	Boolean noting whether the elements in the boundary layerhave free neighbours
	 * 
	 */
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

	
	
	/**
	 * \brief Sets the value of a grid space in the boundary layer, indicating whether the space is part of the boundary layer, or biomass is present
	 * 
	 * Sets the value of a grid space in the boundary layer, indicating whether the space is part of the boundary layer, or biomass is present
	 *  
	 * @param n	The N coordinate of the grid to check whether this square is in the boundary
	 * @param m	The M coordinate of the grid to check whether this square is in the boundary
	 * @param l	The L coordinate of the grid to check whether this square is in the boundary
	 * @return	Integer noting whether or not the square is in the boundary (1 if yes, 0 if not)
	 */
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
						// 2D case  - 
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

	/**
	 * \brief For cyclic boundaries, returns the index of the grid space on the opposite side of the boundary
	 * 
	 * For cyclic boundaries, returns the index of the grid space on the opposite side of the boundary
	 * 
	 * @param val	The integer grid spqce to check
	 * @param limit	The limit of the grid
	 * @return	The integer of the grid square the opposite side of the boundary
	 */
	protected final int cyclicIndex(int val, int limit) {
		return (val<0 ? limit+val : (val>=limit ? val-limit : val));
	}

	/**
	 * \brief Return longest side of this domain
	 * 
	 * Return longest side of this domain
	 * 
	 * @return	Double of the length of the longest side of this domain
	 */
	public double getLongestSize() {
		return Math.max(Math.max(length_X, length_Y), length_Z);
	}

	/**
	 * \brief Return the resolution of this domain
	 * 
	 * Return the resolution of this domain
	 * 
	 * @return	Double value stating the resolution of this domain
	 */
	public double getResolution() {
		return _resolution;
	}

	/**
     * \brief Returns the domain grid
     * 
     * Returns the domain grid
     * 
     * @return SpatialGrid within this domain
     */

	public SpatialGrid getGrid() {
		return _domainGrid;
	}

	/**
	 * \brief Return the name of this domain
	 * 
	 * Return the name of this domain
	 * 
	 * @return Name of this domain
	 */
	public String getName() {
		return domainName;
	}

	/**
	 * \brief Determine if the simulation is recreating a 3D environment
	 * 
	 * Determine if the simulation is recreating a 3D environment
	 * 
	 * @return	A boolean value stating whether or not the environment is 3D
	 */
	public boolean is3D() {
		return _domainGrid.is3D();
	}

	/**
	 * \brief Return the diffusivity grid associated with this domain
	 * 
	 * Return the diffusivity grid associated with this domain
	 * 
	 * @return	SoluteGrid containing diffusivity grid statistics
	 */
	public SoluteGrid getDiffusivity() {
		return _diffusivityGrid;
	}

	/**
	 * \brief Return the boundary layer grid associated with this domain
	 * 
	 * Return the boundary layer associated with this domain
	 * 
	 * @return	SoluteGrid containing boundary between bulk and biofilm
	 */
	public SoluteGrid getBoundaryLayer() 
	{
		return _boundaryLayer;
	}

	/**
	 * \brief Return the biomass grid associated with this domain
	 * 
	 * Return the biomass grid associated with this domain
	 * 
	 * @return	SoluteGrid containing biomass throughout this domain
	 */
	public SoluteGrid getBiomass() {
		return _biomassGrid;
	}

	
	/**
	 * \brief Used in testing to check the top of the boundary layer was calculated correctly
	 * 
	 * Used in testing to check the top of the boundary layer was calculated correctly. KA 210513
	 */
	public void printTopOfBoundaryLayerArray()
	{
		for(int k=1;k<_nK+1;k++)
		{
			for(int j=1;j<_nJ+1;j++)
			{
				System.out.print(_topOfBoundaryLayer[j][k]+" ");
			}
			System.out.println();
		}
	}
	
	/**
	 * \brief Used in testing to view the boundary layer matrix for a set part of the domain
	 * 
	 * Used in testing to view the boundary layer matrix for a set part of the domain. KA 210513
	 */
	public void printBoundaryLayer()
	{
		// Printing the Boundary Layer Grid
	
		for(int k=1;k<_boundaryLayer.getGridSizeK()+1;k++)
		{
			for(int i=1;i<_boundaryLayer.getGridSizeI()+1;i++)
			{
				for(int j=1;j<_boundaryLayer.getGridSizeJ()+1;j++)
				{
					System.out.print(_boundaryLayer.getValueAt(i, j, k)+" ");
				}
				System.out.println();
			}
		}
	}
	
	/**
	 * \brief Used in testing to view the biomass matrix for a set part of the domain
	 * 
	 * Used in testing to view the biomass matrix for a set part of the domain. KA 210513
	 */
	public void printBioMassGrid()
	{
		// Printing the Boundary Layer Grid
	
		for(int k=1;k<_biomassGrid.getGridSizeK()+1;k++)
		{
			for(int i=1;i<_biomassGrid.getGridSizeI()+1;i++)
			{
				for(int j=1;j<_biomassGrid.getGridSizeJ()+1;j++)
				{
					System.out.print(_biomassGrid.getValueAt(i, j, k)+" ");
				}
				System.out.println();
			}
		}
	}
	
	// bvm note 16.12.08: copied this routine following above for the boundary
		// returns true only if there is a free neighbor in the biomass grid
		/*private boolean bflmHasFreeNbh() {
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
		}*/
}
