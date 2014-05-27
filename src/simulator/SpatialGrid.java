/**
 * \package simulator
 * \brief Package of classes that create a simulator object and capture simulation time.
 * 
 * Package of classes that create a simulator object and capture simulation time. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator;

import java.io.Serializable;
import java.util.Arrays;

import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.agent.LocatedAgent;

import utils.ExtraMath;
import utils.MatrixOperations;
import utils.ResultFile;

/**
 * \brief Class defining a spatial grid, i.e. a matrix of double. The grid is padded, 3D grid
 * 
 * Class defining a spatial grid, i.e. a matrix of double. The grid is padded, 3D grid
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class SpatialGrid implements Serializable {

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Name assigned to this spacial grid. Taken from an XML tag in the protocol file
	 */
	public String             gridName;
	
	/**
	 * The unit for all values stored on this grid
	 */
	public String			  gridUnit = "g.L-1";

	/**
	 * The solute grid - a three dimesional array of double values
	 */
	public double[][][]       grid;

	/**
	 * Number of grid voxels in I direction
	 */
	protected int             _nI;
	
	/**
	 * Number of grid voxels in J direction
	 */
	protected int 			  _nJ;
	
	/**
	 * Number of grid voxels in K direction
	 */
	protected int			  _nK;

	/**
	 * Grid resolution = side length of a voxel
	 */
	protected double          _reso;

	/**
	 * Boolean noting whether this grid is 3D (true) or 2D (false)
	 */
	protected boolean         _is3D;

	
	/**
	 * \brief Blank constructor
	 * 
	 * Blank constructor
	 */
	public SpatialGrid() {
	}

	/**
	 * \brief Default constructor for an empty spatial grid
	 * 
	 * Default constructor for an empty spatial grid. Sets the grid resolution and dimensions as provided in the simulation protocol file
	 * 
	 * @param nI	The number of grid locations in the I direction
	 * @param nJ	The number of grid locations in the J direction
	 * @param nK	The number of grid locations in the K direction
	 * @param resolution the grid resolution
	 */
	public SpatialGrid(int nI, int nJ, int nK, double resolution) {

		//sonia:chemostat
		// set the size of the spatial grid to 1,1,1
		if(Simulator.isChemostat){
			_nI = 1;
			_nJ = 1;
			_nK = 1;
		}else{
			_nI = nI;
			_nJ = nJ;
			_nK = nK;
		}
		_reso = resolution;
		// Create a padded grid
		initGrids();

	}

	/**
	 * \brief Default constructor for an empty 2D spatial grid
	 * 
	 * Default constructor for an empty 2D spatial grid. Sets the grid resolution and dimensions as provided in the simulation protocol file
	 * 
	 * @param nI	The number of grid locations in the I direction
	 * @param nJ	The number of grid locations in the J direction
	 * @param resolution the grid resolution
	 */
	public SpatialGrid(int nI, int nJ, double resolution) {
		//sonia:chemostat
		// Rob 17/8/2011: this doesn't seem to be called by anything

		if(Simulator.isChemostat){
			_nI = 1;
			_nJ = 1;
		}else{    // 2D biofilm
			_nI = nI;
			_nJ = nJ;
		}
		_nK = 1;
		_reso = resolution;
		_is3D = false;
		// Create a padded grid
		initGrids();
		
	}

	/**
	 * \brief Creates the solute grid at the required size
	 * 
	 * Creates the solute grid at the required size. If this is a chemostat, this will simple be a 1x1x1, if not there are further
	 * checks to determine whether we are simulating 3D or not
	 */
	protected void initGrids() 
	{
		if(Simulator.isChemostat)
		{
			//sonia:chemostat
			//set the size of the grid to be 1,1,1 and not a padded grid
			
			_is3D = false;
			grid = new double [1][1][1];

		}
		else
		{
			// Obviously if we create only one cell in the Z dimension, the grid is 2D
			_is3D = _nK > 1;

			grid = new double[_nI+2][_nJ+2][_nK+2];
		}
	
		// KA 28/3/13 - this was turned off when I started recommenting iDynoMiCS - no reason given here but may need checking
		// At the creation the table is automatically filled by 0
		//setAllValueAt(0);

	}

	/**
	 * \brief Determine if a given discrete position is valid or outside the grid (Padding excluded)
	 * 
	 * Determine if a given discrete position is valid or outside the grid (Padding excluded)
	 * 
	 * @param dC	DiscreteVector to validate
	 * @return Boolean stating whether this location is valid (true) or outside the grid
	 */
	public boolean isValid(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i>=0)&(dC.i<_nI);
		out &= (dC.j>=0)&(dC.j<_nJ);
		out &= (dC.k>=0)&(dC.k<_nK);
		return out;
	}

	/**
	 * \brief Determine if a given discrete position is valid or outside the grid (Padding included)
	 * 
	 * Determine if a given discrete position is valid or outside the grid
	 * 
	 * @param dC	DiscreteVector to validate
	 * @return Boolean stating whether this location is valid (true) or outside the grid
	 */
	public boolean isValidOrPadded(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i>=-1)&(dC.i<=_nI);
		out &= (dC.j>=-1)&(dC.j<=_nJ);
		out &= (dC.k>=-1)&(dC.k<=_nK);
		return out;
	}

	/**
	 * \brief Determine if a given voxel coordinate is valid or outside the grid
	 * 
	 * Determine if a given voxel coordinate is valid or outside the grid
	 * 
	 * @param i	I Coordinate of the grid location
	 * @param j	J Coordinate of the grid location
	 * @param k K Coordinate of the grid location
	 * @return Boolean stating whether this location is valid (true) or outside the grid
	 */
	public boolean isValidorPadded(int i, int j, int k) {
		boolean out = true;
		out &= (i>=0)&(i<=_nI);
		out &= (j>=0)&(j<=_nJ);
		out &= (k>=0)&(k<=_nK);
		return out;
	}

	/**
	 * \brief Determine if a given continuous location is valid or outside the grid
	 * 
	 * Determine if a given continuous location is valid or outside the grid
	 * 
	 * @param cc	ContinuousVector to validate
	 * @return Boolean stating whether this location is valid (true) or outside the grid
	 */
	public boolean isValid(ContinuousVector cc) {
		return isValid(getDiscreteCoordinates(cc));
	}

	/**
	 * \brief Test if a location is valid or is outside the grid. Note padding is included as valid
	 * 
	 * Test if a location is valid or is outside the grid. Note padding is included as valid
	 * 
	 * @param cc	ContinuousVector containing the location to be tested
	 * @return	Boolean stating whether this location is valid (true) or outside the grid
	 */
	public boolean isValidOrPadded(ContinuousVector cc) {
		return isValidOrPadded(getDiscreteCoordinates(cc));
	}

	/**
	 * \brief Transform a location, expressed as a continuous vector into a discrete position on the basis of the resolution of the grid
	 * 
	 * Transform a location, expressed as a continuous vector into a discrete position on the basis of the resolution of the grid
	 * 
	 * @param cC	ContinuousVector to be transformed
	 * @return	DiscreteVector created from this continuous location
	 */
	public DiscreteVector getDiscreteCoordinates(ContinuousVector cC) {
		int i = (int) Math.floor(cC.x/_reso);
		int j = (int) Math.floor(cC.y/_reso);
		int k = (int) Math.floor(cC.z/_reso);

		return new DiscreteVector(i, j, k);
	}

	/**
	 * \brief Transform a position, expressed as a discrete vector into a continuous location on the basis of the resolution of the grid
	 * 
	 * Transform a position, expressed as a discrete vector into a continuous location on the basis of the resolution of the grid
	 * 
	 * @param dC	DiscreteVector to be transformed
	 * @return	ContinuousVector created from this discrete position
	 */
	public ContinuousVector getContinuousCoordinates(DiscreteVector dC) {
		ContinuousVector cc = new ContinuousVector(dC, _reso);
		return cc;
	}

	/**
	 * \brief Convert a set of voxel coordinates into a continuous vector
	 * 
	 * Convert a set of voxel coordinates into a continuous vector
	 * 
	 * @param i	I Coordinate of the grid location
	 * @param j	J Coordinate of the grid location
	 * @param k K Coordinate of the grid location
	 * @return	ContinuousVector created from these voxel coordinates
	 */
	public ContinuousVector getContinuousCoordinates(int i, int j, int k) {
		ContinuousVector cc = new ContinuousVector((i+.5)*_reso, (j+.5)*_reso, (k+.5)*_reso);
		return cc;
	}

	/**
	 * \brief Return the maximum value on this grid (padding included)
	 * 
	 * Return the maximum value on this grid (padding included)
	 * 
	 * @return Maximum value of the grid
	 */
	public double getMax() {
		return MatrixOperations.max(grid);
	}

	/**
	 * \brief Return the average value on this grid (padding included)
	 * 
	 * Return the average value on this grid (padding included)
	 * 
	 * @return Average value of the grid
	 */
	public double getAverage() {
		return MatrixOperations.computeSum(grid)/(_nI)/(_nJ)/(_nK);
	}

	/**
	 * \brief Return the sum of this grid
	 * 
	 * Return the sum of this grid
	 * 
	 * @return	The sum of the values in this spatial grid
	 */
	public double getSum() {
		return MatrixOperations.computeSumP(grid);
	}

	/**
	 * \brief Return the minimum value on this grid (padding included)
	 * 
	 * Return the minimum value on this grid (padding included)
	 * 
	 * @return Minimum value of the grid
	 */
	public double getMin() {
		return MatrixOperations.min(grid);
	}

	/**
	 * \brief For a given location, calculate the 2nd spatial derivative according to X
	 * 
	 * For a given location, calculate the 2nd spatial derivative according to X
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 2nd spatial derivative according X
	 */
	public double diff2X(int i, int j, int k) {
		double d2x = (grid[i+1][j][k] + grid[i-1][j][k] - 2 * grid[i][j][k]) / ExtraMath.sq(_reso);
		return (Double.isNaN(d2x) || Double.isInfinite(d2x)) ? 0 : d2x;
	}

	/**
	 * \brief For a given location, expressed as a discrete vector, calculate the 2nd spatial derivative according to X
	 * 
	 * For a given location, expressed as a discrete vector, calculate the 2nd spatial derivative according to X
	 * 
	 * @param dV	DiscreteVector containing the position of a grid location
	 * @return 2nd spatial derivative according X
	 */
	public double diff2X(DiscreteVector dV) {
		double d2x = (grid[dV.i+1][dV.j][dV.k]
                        + grid[dV.i-1][dV.j][dV.k]
                        - 2 * grid[dV.i][dV.j][dV.k]
                     ) / ExtraMath.sq(_reso);
		return (Double.isNaN(d2x) || Double.isInfinite(d2x)) ? 0 : d2x;
	}

	/**
	 * \brief For a given location, calculate the 1st spatial derivative according to X
	 * 
	 * For a given location, calculate the 1st spatial derivative according to X
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 1st spatial derivative according X
	 */
	public double diffX(int i, int j, int k) {
		double dx = (grid[i+1][j][k] - grid[i-1][j][k]) / (2*_reso);
		return (Double.isNaN(dx) || Double.isInfinite(dx)) ? 0 : dx;
	}

	/**
	 * \brief For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to X
	 * 
	 * For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to X
	 * 
	 * @param dV	DiscreteVector containing the position of a grid location
	 * @return 1st spatial derivative according X
	 */
	public double diffX(DiscreteVector dV) {
		double dx = (grid[dV.i+1][dV.j][dV.k] - grid[dV.i-1][dV.j][dV.k]) / (2*_reso);
		return (Double.isNaN(dx) || Double.isInfinite(dx)) ? 0 : dx;
	}

	/**
	 * \brief For a given location, calculate the 2nd spatial derivative according to Y
	 * 
	 * For a given location, calculate the 2nd spatial derivative according to Y
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 2nd spatial derivative according Y
	 */
	public double diff2Y(int i, int j, int k) {
		double d2y =  (grid[i][j+1][k] + grid[i][j-1][k] - 2 * grid[i][j][k]) / ExtraMath.sq(_reso);
		return (Double.isNaN(d2y) || Double.isInfinite(d2y)) ? 0 : d2y;
	}

	/**
	 * \brief For a given location, expressed as a discrete vector, calculate the 2nd spatial derivative according to Y
	 * 
	 * For a given location, expressed as a discrete vector, calculate the 2nd spatial derivative according to Y
	 * 
	 * @param dV	DiscreteVector containing the position of a grid location
	 * @return 2nd spatial derivative according Y
	 */
	public double diff2Y(DiscreteVector dV) {
		double d2y =  (grid[dV.i][dV.j+1][dV.k]
                        + grid[dV.i][dV.j-1][dV.k]
                        - 2 * grid[dV.i][dV.j][dV.k]
                      ) / ExtraMath.sq(_reso);
		return (Double.isNaN(d2y) || Double.isInfinite(d2y)) ? 0 : d2y;
	}

	/**
	 * \brief For a given location, calculate the 1st spatial derivative according to Y
	 * 
	 * For a given location, calculate the 1st spatial derivative according to Y
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 1st spatial derivative according Y
	 */
	public double diffY(int i, int j, int k) {
		double dy =  (grid[i][j+1][k] - grid[i][j-1][k]) / (2*_reso);
		return (Double.isNaN(dy) || Double.isInfinite(dy)) ? 0 : dy;
	}
	
	/**
	 * \brief For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to Y
	 * 
	 * For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to Y
	 * 
	 * @param dV	DiscreteVector containing the position of a grid location
	 * @return 1st spatial derivative according Y
	 */
	public double diffY(DiscreteVector dV) {
		double dy =  (grid[dV.i][dV.j+1][dV.k] - grid[dV.i][dV.j-1][dV.k]) / (2*_reso);
		return (Double.isNaN(dy) || Double.isInfinite(dy)) ? 0 : dy;
	}
	

	/**
	 * \brief For a given location, calculate the 2nd spatial derivative according to Z
	 * 
	 * For a given location, calculate the 2nd spatial derivative according to Z
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 2nd spatial derivative according Z
	 */
	public double diff2Z(int i, int j, int k) {
		double d2z = (grid[i][j][k+1] + grid[i][j][k-1] - 2 * grid[i][j][k]) / (_reso*_reso);
		return (Double.isNaN(d2z) || Double.isInfinite(d2z)) ? 0 : d2z;
	}

	/**
	 * \brief For a given location, calculate the 1st spatial derivative according to Z
	 * 
	 * For a given location, calculate the 1st spatial derivative according to Z
	 * 
	 * @param i	I position on the spatial grid
	 * @param j	J position on the spatial grid
	 * @param k	K position on the spatial grid
	 * @return 1st spatial derivative according Z
	 */
	public double diffZ(int i, int j, int k) {
		double dz =  (grid[i][j][k+1] - grid[i][j][k-1]) / (2*_reso);
		return (Double.isNaN(dz) || Double.isInfinite(dz)) ? 0 : dz;
	}
	
	/**
	 * \brief For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to Z
	 * 
	 * For a given location, expressed as a discrete vector, calculate the 1st spatial derivative according to Z
	 * 
	 * @param dV	DiscreteVector containing the position of a grid location
	 * @return 1st spatial derivative according Z
	 */
	public double diffZ(DiscreteVector dV) {
		double dz =  (grid[dV.i][dV.j][dV.k+1] - grid[dV.i][dV.j][dV.k-1]) / (2*_reso);
		return (Double.isNaN(dz) || Double.isInfinite(dz)) ? 0 : dz;
	}
	/**
	 * \brief Computes the average concentration seen in a sphere (or cube) centered around a given point
	 * 
	 * Computes the average concentration seen in a sphere (or cube) centered around a given point
	 * 
	 * @param cC	ContinuousVector containing the point to use as the centre of this search
	 * @param extReso	Resolution to use in this search
	 * @return	Average grid value seen around this point
	 */
	public double getValueAround(ContinuousVector cC, double extReso) {

		if (extReso<=_reso) {
			// The asked value is
			return getValueAt(cC);
		} else {
			// TODO
			return getValueAt(cC);
		}
	}


	

	/**
	 * \brief Returns the average grid value seen in a sphere (or cube) located around a given agent
	 * 
	 * Returns the average grid value seen in a sphere (or cube) located around a given agent
	 * 
	 * @param aLocAgent	The agent to use as the centre of the search
	 * @return	Average grid value seen around this point
	 */
	public double getValueAround(LocatedAgent aLocAgent) {
		return getValueAround(aLocAgent.getLocation(), aLocAgent.getRadius(true));
	}


	/**
	 * \brief Return the value stored at the location given by the stated continuous vector
	 * 
	 * Return the value stored at the location given by the stated continuous vector
	 * 
	 * @param cC	ContinuousVector containing the grid location to return
	 * @return	Double value stored at that grid location
	 */
	public double getValueAt(ContinuousVector cC) {
		DiscreteVector dc = getDiscreteCoordinates(cC);
		//sonia:chemostat

		if(Simulator.isChemostat){
			if (isValid(dc)) {
				return grid[dc.i][dc.j][dc.k];
			} else {
				return Double.NaN;
			}
		}else{

			if (isValid(dc)) 
			{
				return grid[dc.i+1][dc.j+1][dc.k+1];
				
			} else {
				return Double.NaN;
			}
		}
	}

	/**
	 * \brief Returns a vector of the first spatial derivatives in x, y & z 
	 * 
	 * Returns a vector of the first spatial derivatives in x, y & z (nabla cC - see http://en.wikipedia.org/wiki/Del). Does this by 
	 * first converting the ContinuousVector to a DiscreteVector and then estimating then gradient using the Mean Value Theorem 
	 * (http://en.wikipedia.org/wiki/Mean_value_theorem).
	 * 
	 * @param cC	ContinuousVector position used to calculate the gradient
	 * @return	Vector of spatial derivatives in X,Y,Z
	 */
	public ContinuousVector getGradient(ContinuousVector cC) {
		DiscreteVector dV = new DiscreteVector(cC,_reso);
		return new ContinuousVector(diffX(dV), diffY(dV), diffZ(dV));
	}

	/**
	 * \brief Returns a vector of the first spatial derivatives in x and y, for 2D simulations 
	 * 
	 * Returns a vector of the first spatial derivatives in x and y (nabla cC - see http://en.wikipedia.org/wiki/Del). Does this by 
	 * first converting the ContinuousVector to a DiscreteVector and then estimating then gradient using the Mean Value Theorem 
	 * (http://en.wikipedia.org/wiki/Mean_value_theorem).
	 * 
	 * @param cC	ContinuousVector position used to calculate the gradient
	 * @return	Vector of spatial derivatives in X and Y
	 */
	public ContinuousVector getGradient2D(ContinuousVector cC) 
	{
		DiscreteVector dV = new DiscreteVector(cC,_reso);
		return new ContinuousVector(diffX(dV), diffY(dV), diffY(dV));
	}

	//public ContinuousVector getGradient(DiscreteVector dV) {
	//	return new ContinuousVector(diffX(dV), diffY(dV), diffZ(dV));
	//}

	/**
	 * \brief Return the value on the padded grid at a given position (the coordinates are NOT corrected)
	 * 
	 * Return the value on the padded grid at a given position (the coordinates are NOT corrected)
	 * 
	 * @param dV	DiscreteVector containing the location of the grid whos value should be returned
	 * @return The double value at that location
	 */
	public double getValueAt(DiscreteVector dV) {
		//sonia:chemostat

		if(Simulator.isChemostat){
			if (isValid(dV)) return grid[dV.i][dV.j][dV.k];
			else return Double.NaN;	
		}else
		{
				if (isValid(dV)) return grid[dV.i+1][dV.j+1][dV.k+1];
				else return Double.NaN;
			
			
		}
	}

	/**
	 * \brief Return the value on the padded grid at a given position (the coordinates are NOT corrected)
	 * 
	 * Return the value on the padded grid at a given position (the coordinates are NOT corrected)
	 * 
	 * @param i	I Coordinate of the grid location to set
	 * @param j	J Coordinate of the grid location to set
	 * @param k K Coordinate of the grid location to set
	 * @return The double value at that location
	 */
	public double getValueAt(int i, int j, int k) 
	{
		
			if (isValidorPadded(i, j, k)) return grid[i][j][k];
			else return Double.NaN;
		
		
	}

	/**
	 * \brief Set a grid location, expressed as a ContinuousVector, to a specified value. The coordinates are corrected
	 * 
	 * Set a grid location, expressed as a ContinuousVector, to a specified value. The coordinates are corrected
	 * 
	 * @param value	Value to set the specified location to
	 * @param cC	Continuous vector stating the location of the grid to be set to the given value
	 */
	public void setValueAt(double value, ContinuousVector cC) {
		DiscreteVector dC = getDiscreteCoordinates(cC);
		setValueAt(value, dC);
	}

	/**
	 * \brief Set a grid location, expressed as a DiscreteVector, to a specified value. The coordinates are corrected
	 * 
	 * Set a grid location, expressed as a DiscreteVector, to a specified value. The coordinates are corrected
	 * 
	 * @param value	Value to set the specified location to
	 * @param dC	Discrete vector stating the location of the grid to be set to the given value
	 */
	public void setValueAt(double value, DiscreteVector dC) 
	{
		//sonia:chemostat
		if(Simulator.isChemostat){
			grid[dC.i][dC.j][dC.k] = value;
		}else
		{
			grid[dC.i+1][dC.j+1][dC.k+1] = value;
		}
	}

	/**
	 * \brief Set a grid location to a specified value. Note the coordinates are NOT corrected
	 * 
	 * Set a grid location to a specified value. Note the coordinates are NOT corrected
	 * 
	 * @param value	Value to set the grid location to
	 * @param i	I Coordinate of the grid location to set
	 * @param j	J Coordinate of the grid location to set
	 * @param k K Coordinate of the grid location to set
	 */
	public void setValueAt(double value, int i, int j, int k) 
	{
		grid[i][j][k] = value;
	}

	/**
	 * \brief Add a value to that contained at the given discrete coordinates of this grid. Coordinates are corrected for padding
	 * 
	 * Add a value to that contained at the given discrete coordinates of this grid. Coordinates are corrected for padding
	 * 
	 * @param value	Value to add to the specified grid location
	 * @param cC	Continuous vector expressing the location of the grid to be increased
	 */
	public void addValueAt(double value, ContinuousVector cC) 
	{
		addValueAt(value, getDiscreteCoordinates(cC));
	}

	/**
	 * \brief Add a value to that contained at the given discrete coordinates of this grid. Coordinates are corrected for padding
	 * 
	 * Add a value to that contained at the given discrete coordinates of this grid. Coordinates are corrected for padding
	 * 
	 * @param value	Value to add to the specified grid location
	 * @param dC	Discrete vector expressing the location of the grid to be increased
	 */
	public void addValueAt(double value, DiscreteVector dC) {
		//sonia:chemostat
		if(Simulator.isChemostat){
			grid[dC.i][dC.j][dC.k] += value;
		}else
		{
			grid[dC.i+1][dC.j+1][dC.k+1] += value;
			
		}
	}

	/**
	 * \brief Add a value to all locations on this grid (including the padding)
	 * 
	 * Add a value to all locations on this grid (including the padding)
	 * 
	 * @param value	Value to be added to the contents of all grid voxels
	 */
	public void addAllValues(double value) {
		//sonia:chemostat

		if(Simulator.isChemostat){
			for (int i = 0; i<_nI; i++) {
				for (int j = 0; j<_nJ; j++) {
					for (int k = 0; k<_nK; k++) {
						grid[i][j][k] += value;
					}
				}
			}

		}else{
			for (int i = 0; i<_nI+2; i++) {
				for (int j = 0; j<_nJ+2; j++) {
					for (int k = 0; k<_nK+2; k++) {
						grid[i][j][k] += value;
					}
				}
			}
		}
	}

	/**
	 * \brief Checks a value at a given location and sets it to zero if the value is negative
	 * 
	 * Checks a value at a given location and sets it to zero if the value is negative
	 * 
	 * @param i	Voxel coordinate in I direction
	 * @param j	Voxel coordinate in J direction
	 * @param k	Voxel coordinate in K direction
	 */
	public void truncateValueAt(int i, int j, int k) {
		grid[i][j][k] = (grid[i][j][k]<0 ? 0 : grid[i][j][k]);
	}

	/**
	 * \brief Set all meshes of a grid with the same value (including the padding - if not a chemostat run)
	 * 
	 * Set all meshes of a grid with the same value (including the padding if not a chemostat simulation)
	 * 
	 * @param value	Value at which to set all the elements of the grid
	 */
	public void setAllValueAt(double value) 
	{
		//sonia:chemostat 
		//in this case we have no padding

		if(Simulator.isChemostat){
			Arrays.fill(grid[0][0],value);

		}else{
			for (int i = 0; i<_nI+2; i++) {
				for (int j = 0; j<_nJ+2; j++) {
					Arrays.fill(grid[i][j], value);
				}
			}
		}
	}


	/**
	 * \brief Return the number of voxels in the X direction
	 * 
	 * Return the number of voxels in the X direction (ignoring the padding)
	 * 
	 * @return Number of real voxels along X 
	 */
	public int getGridSizeI() 
	{
		return _nI;
	}
	
	/**
	 * \brief Return the number of voxels in the Y direction
	 * 
	 * Return the number of voxels in the Y direction (ignoring the padding)
	 * 
	 * @return Number of real voxels along Y 
	 */
	public int getGridSizeJ() 
	{
		return _nJ;
	}

	/**
	 * \brief Return the number of voxels in the Z direction
	 * 
	 * Return the number of voxels in the Z direction (ignoring the padding)
	 * 
	 * @return Number of real voxels along Z 
	 */
	public int getGridSizeK() 
	{
		return _nK;
	}

	/**
	 * \brief Return the number of voxels along a given direction (axeCode)
	 * 
	 * Return the number of voxels along a given direction (axeCode):  1-X, 2-Y,3-Z
	 * 
	 * @param axeCode Integer noting the direction to query
	 * @return : the number of voxels along a direction including padding bands
	 */
	public int getGridTotalSize(int axeCode) 
	{
		//sonia:chemostat
		if(Simulator.isChemostat){
			switch (axeCode) {
			case 1:
				return _nI;
			case 2:
				return _nJ;
			case 3:
				return _nK;
			default:
				return 0;
			}
		}else{
			switch (axeCode) {
			case 1:
				return _nI+2;
			case 2:
				return _nJ+2;
			case 3:
				return _nK+2;
			default:
				return 0;
			}
		}
	}

	/**
	 * \brief Returns the length (in distance unit) along a given direction
	 * 
	 * Returns the length (in distance unit) along a given direction. The direction is stated as an integer: 1 for X, 2 for Y, 3 for Z
	 * 
	 * @param axeCode	The direction of which the length is required: 1-X, 2-Y, 3-Z
	 * @return Double value stating the length (in distance unit) along a direction ignoring padding bands
	 */ 
	public double getGridLength(int axeCode) 
	{
		switch (axeCode) {
		case 1:
			return _nI*_reso;
		case 2:
			return _nJ*_reso;
		case 3:
			return _nK*_reso;
		default:
			return (double) 0;
		}
	}

	/**
	 * \brief Return the volume of one voxel of the spatial grid
	 * 
	 * Return the volume of one voxel of the spatial grid
	 * 
	 * @return	Double value stating the volume of one voxel of the spatial grid
	 */
	public double getVoxelVolume() {
		return ExtraMath.cube(_reso);
	}

	/**
	 * \brief Return the whole grid including the padding
	 * 
	 * Return the whole grid including the padding
	 * 
	 * @return the spatial grid
	 */
	public double[][][] getGrid() {
		return grid;
	}

	/**
	 * \brief Return a clone of this spatial grid
	 * 
	 * Return a clone of this spatial grid
	 * 
	 * @return	A clone of this spatial grid
	 */
	public double[][][] getCloneGrid(){
		return grid.clone();
	}

	/**
	 * \brief Returns the resolution of this spatial grid
	 * 
	 * Returns the resolution of this spatial grid
	 * 
	 * @return	Double value stating the resolution (in micrometers) of this grid
	 */
	public double getResolution() 
	{
		return _reso;
	}

	/**
	 * \brief Determine if this spatial grid is 3D or 2D
	 * 
	 * Determine if this spatial grid is 3D or 2D
	 * 
	 * @return	Boolean noting whether this grid is 3D (true) or 2D (false)
	 */
	public boolean is3D() {
		//sonia:chemostat
		if(Simulator.isChemostat){
			return false;
		}else{
			return _is3D;
		}
	}

	/**
	 * \brief Set the values of this spatial grid to those contained in the supplied grid
	 * 
	 * @param u	Matrix of values which to set the spatial grid to
	 */
	public void setGrid(double[][][] u) {
		utils.MatrixOperations.copyValuesTo(grid,u);
	}

	/**
	 * \brief Write the contents of this grid to the XML results files
	 * 
	 * Write the contents of this grid to the XML results files. This shows the level of solute in each of the grid spaces
	 * 
	 * @param bufferState	The output buffer writing the env_state file for this iteration
	 * @param bufferSummary	The output buffer writing the env_sum file for this iteration
	 * @throws Exception	Exception thrown if these buffers cannot be opened for writing to
	 */
	public void writeReport(ResultFile bufferState, ResultFile bufferSummary) throws Exception {

		// Edit the markup for the solute grid
		StringBuffer value = new StringBuffer();
		value.append("<solute name=\"").append(gridName);
		value.append("\" unit=\"").append(gridUnit);
		value.append("\" resolution=\"").append(_reso);
		value.append("\" nI=\"").append(_nI);
		value.append("\" nJ=\"").append(_nJ);
		value.append("\" nK=\"").append(_nK);
		value.append("\">\n");

		// Write the markup in the file
		bufferState.write(value.toString());

		// Rob 3/3/11: Changed to fix bug in envState output files
		// and improve code readability (plus efficiency... possibly).
		// Note that for a chemostat, i=j=k=1 and that in 2D k=1. The
		// main fix here however, is that grid is a double and not an
		// array, as previously coded (this reduces the amount of 
		// storage space taken by envState files by about a 2 thirds!)

		//sonia:chemostat
		if(Simulator.isChemostat){
			// bufferState.write(Arrays.toString(grid[i][j]));
			bufferState.write(Double.toString(grid[0][0][0]));
			bufferState.write(";\n");

		}else{

			// KA 06062013 - turned off the printing of the padding. Will need to ensure this is clear from v1.2
			// Fill the mark-up
			if (_nK==1) {
				// We have a 2D grid
				
				
				for (int i = 1; i<_nI+1; i++) {
					for (int j = 1; j<_nJ+1; j++) {
						//bufferState.write(Arrays.toString(grid[i][j]));
						bufferState.write(Double.toString(grid[i][j][1]));
						bufferState.write(";\n");
					}
				}
			} else {
				// We have a 3D grid 
				for (int i = 1; i<_nI+1; i++) {
					for (int j = 1; j<_nJ+1; j++) {
						for (int k=1; k<_nK+1; k++) {
							// bufferState.write(Arrays.toString(grid[i][j]));
							bufferState.write(Double.toString(grid[i][j][k]));
							bufferState.write(";\n");
						}
					}
				}
			}
		}

		// Close the mark-up
		bufferState.write("\n</solute>\n");
		
		
		

	}

}
