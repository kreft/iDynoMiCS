/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 * ________________________________________________________
 * Class defining a spatial grid, i.e. a matrix of double.
 * The grid is padded, 3D grid
 */

/**
 * @since june 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
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
import utils.XMLParser;

public class SpatialGrid implements Serializable {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	public String             gridName;
	public String			  gridUnit = "g.L-1";

	// the grid
	public double[][][]       grid;

	// size of original grid
	protected int             _nI, _nJ, _nK;

	// grid resolution = side length of a voxel
	protected double          _reso;

	protected boolean         _is3D;

	/* ______________________________ CONSTRUCTOR ___________________________ */

	public SpatialGrid() {
	}

	/**
	 * Default constructor for an empty array
	 * 
	 * @param resolution the grid resolution
	 * @param nI,nJ,nK the original grid size
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

	public SpatialGrid(XMLParser cdRoot) {
		//sonia:chemostat
		// Rob (17/8/2011): this doesn't seem to be called by anything 

		if(Simulator.isChemostat){
			_nI = 1;
			_nJ = 1;
			_nK = 1;
		}else{

			// First test number of dimensions (2 or 3)
			_is3D = (cdRoot.getAttributeDbl("nDim")==3);

			_nI = (int) cdRoot.getAttributeDbl("nI");
			_nJ = (int) cdRoot.getAttributeDbl("nJ");
			if (_is3D) {
				_nK = (int) cdRoot.getAttributeDbl("nK");
			} else {
				_nK = 1;
			}
		}
		_reso = cdRoot.getParamLength("resolution");

		// Create a padded grid
		//sonia:chemostat
		//I've changed initGrids() to avoid padding
		initGrids();

	}

	/**
	 * Empty 2D-grid buidler
	 * 
	 * @param nI
	 * @param nJ
	 * @param res
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
	 * Create the grid and fill the value to 0
	 */
	protected void initGrids() {
		//sonia:chemostat
		//set the size of the grid to be 1,1,1 and not a padded grid
		if(Simulator.isChemostat){
			_is3D = false;
			grid = new double [1][1][1];

		}else{
			// Obviously if we create only one cell in the Z dimension, the grid is
			// 2D
			_is3D = !(_nK==1);

			grid = new double[_nI+2][_nJ+2][_nK+2];
		}
		// At the creation the table is automatically filled by 0
		// setAllValueAt(0);

	}

	/* _______________________ TOOLS __________________________________ */

	/**
	 * Test if a position is defined or outside the grid, here padding bands are
	 * not considered as valid coordinates
	 * 
	 * @param cc
	 * @return
	 */
	public boolean isValid(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i>=0)&(dC.i<_nI);
		out &= (dC.j>=0)&(dC.j<_nJ);
		out &= (dC.k>=0)&(dC.k<_nK);
		return out;
	}

	/**
	 * Test if a position is defined or outside the grid, padding bands are now
	 * valid coordinates
	 * 
	 * @param cc
	 * @return
	 */
	public boolean isValidOrPadded(DiscreteVector dC) {
		boolean out = true;
		out &= (dC.i>=-1)&(dC.i<=_nI);
		out &= (dC.j>=-1)&(dC.j<=_nJ);
		out &= (dC.k>=-1)&(dC.k<=_nK);
		return out;
	}

	public boolean isValidorPadded(int i, int j, int k) {
		boolean out = true;
		out &= (i>=0)&(i<=_nI);
		out &= (j>=0)&(j<=_nJ);
		out &= (k>=0)&(k<=_nK);
		return out;
	}

	/**
	 * Test if a location is defined or outside the grid, padding bands are not
	 * considered as valid coordinates
	 * 
	 * @param cc
	 * @return
	 */
	public boolean isValid(ContinuousVector cc) {
		return isValid(getDiscreteCoordinates(cc));
	}

	/**
	 * Test if a location is defined or outside the grid, padding bands are now
	 * considered as valid coordinates
	 * 
	 * @param cc
	 * @return
	 */
	public boolean isValidOrPadded(ContinuousVector cc) {
		return isValidOrPadded(getDiscreteCoordinates(cc));
	}

	public void readGrid(SpatialGrid aSpG) {
		for (int i = 0; i<_nI+2; i++) {
			for (int j = 0; j<_nJ+2; j++) {
				for (int k = 0; k<_nK+2; k++) {
					grid[i][j][k] = aSpG.grid[i][j][k];
				}
			}
		}
	}

	/* ____________________ USEFULL TOOLS ______________________________ */

	/**
	 * Transform a location (continuous) in a position (discrete) on the basis
	 * of the resolution of the grid
	 */
	public DiscreteVector getDiscreteCoordinates(ContinuousVector cC) {
		int i = (int) Math.floor(cC.x/_reso);
		int j = (int) Math.floor(cC.y/_reso);
		int k = (int) Math.floor(cC.z/_reso);

		return new DiscreteVector(i, j, k);
	}

	/**
	 * Transform a position (discrete) in a location (continuous) on the basis
	 * of the resolution of the grid
	 */
	public ContinuousVector getContinuousCoordinates(DiscreteVector dC) {
		ContinuousVector cc = new ContinuousVector(dC, _reso);
		return cc;
	}

	public ContinuousVector getContinuousCoordinates(int i, int j, int k) {
		ContinuousVector cc = new ContinuousVector((i+.5)*_reso, (j+.5)*_reso, (k+.5)*_reso);
		return cc;
	}

	/**
	 * @return the maximal value on a grid (padding bande included)
	 */
	public double getMax() {
		return MatrixOperations.max(grid);
	}

	/**
	 * @return the average value on the grid padding is excluded
	 */
	public double getAverage() {
		return MatrixOperations.computeSum(grid)/(_nI)/(_nJ)/(_nK);
	}

	/*//sonia:chemostat 19.02.2010
	public double getAverageChemo() {
		return MatrixOperations.computeSumChemo(grid)/(_nI)/(_nJ)/(_nK);
	}*/


	public double getSum() {
		return MatrixOperations.computeSum(grid);
	}

	/**
	 * @return the value divided by the surface (in L-2)
	 */
	public double getAreaConc() {
		return MatrixOperations.computeSum(grid)/(_nJ*_reso*_nK*_reso);
	}

	/**
	 * @return the minimum value of the grid (padding bande included)
	 */
	public double getMin() {
		return MatrixOperations.min(grid);
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 2nd spatial derivative according X
	 */
	public double diff2X(int i, int j, int k) {
		double value = (grid[i+1][j][k]+grid[i-1][j][k]-2*grid[i][j][k])/ExtraMath.sq(_reso);		
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	public double diff2X(DiscreteVector dV) {
		double value = (grid[dV.i+1][dV.j][dV.k]+grid[dV.i-1][dV.j][dV.k]-2*grid[dV.i][dV.j][dV.k])/ExtraMath.sq(_reso);		
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 1st spatial derivative according X axe
	 */	
	public double diffX(int i, int j, int k) {
		double value = (grid[i+1][j][k]-grid[i-1][j][k])/(2*_reso);		
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	public double diffX(DiscreteVector dV) {
		double value = (grid[dV.i+1][dV.j][dV.k]-grid[dV.i-1][dV.j][dV.k])/(2*_reso);		
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 2nd spatial derivative according Y axe
	 */
	public double diff2Y(int i, int j, int k) {
		double value =  (grid[i][j+1][k]+grid[i][j-1][k]-2*grid[i][j][k])/ExtraMath.sq(_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	public double diff2Y(DiscreteVector dV) {
		double value =  (grid[dV.i][dV.j+1][dV.k]+grid[dV.i][dV.j-1][dV.k]-2*grid[dV.i][dV.j][dV.k])/ExtraMath.sq(_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 1st spatial derivative according Y axe
	 */
	public double diffY(int i, int j, int k) {
		double value =  (grid[i][j+1][k]-grid[i][j-1][k])/(2*_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}
	public double diffY(DiscreteVector dV) {
		double value =  (grid[dV.i][dV.j+1][dV.k]-grid[dV.i][dV.j-1][dV.k])/(2*_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}
	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 2nd spatial derivative according Z axe
	 */
	public double diff2Z(int i, int j, int k) {
		double value =  (grid[i][j][k+1]+grid[i][j][k-1]-2*grid[i][j][k])/(_reso*_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return 1st spatial derivative according Z axe
	 */
	public double diffZ(int i, int j, int k) {
		double value =  (grid[i][j][k+1]-grid[i][j][k-1])/(2*_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}
	public double diffZ(DiscreteVector dV) {
		double value =  (grid[dV.i][dV.j][dV.k+1]-grid[dV.i][dV.j][dV.k-1])/(2*_reso);
		return ((Double.isNaN(value)|Double.isInfinite(value))?0:value);
	}
	/**
	 * Computes the average concentration seen in a sphere (or cube) centered in
	 * cc
	 * @param cc
	 * @param res
	 * @return
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

	// TODO
	public void setValueAround(double value, ContinuousVector cC, double extReso) {
		setValueAt(value, cC);
	}

	// TODO
	public void addValueAround(double value, ContinuousVector cC, double extReso) {
		addValueAt(value, cC);
	}

	/**
	 * Send the average value saw by a located agent on this grid
	 * 
	 * @param aLocAgent
	 * @return
	 */
	public double getValueAround(LocatedAgent aLocAgent) {
		return getValueAround(aLocAgent.getLocation(), aLocAgent.getRadius(true));
	}

	/**
	 * Set the average value saw by a located agent on this grid
	 * 
	 * @param aLocAgent
	 * @return
	 */
	public void setValueAround(double value, LocatedAgent aLocAgent) {
		setValueAround(value, aLocAgent.getLocation(), aLocAgent.getRadius(true));
	}

	/**
	 * Add the average value saw by a located agent on this grid
	 * 
	 * @param aLocAgent
	 * @return
	 */
	public void addValueAround(double value, LocatedAgent aLocAgent) {
		addValueAround(value, aLocAgent.getLocation(), aLocAgent.getRadius(true));
	}

	/* ______________________ GET & SET _____________________________________ */

	/**
	 * Send the value on the padded grid at a given location
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

			if (isValid(dc)) {
				return grid[dc.i+1][dc.j+1][dc.k+1];
			} else {
				return Double.NaN;
			}
		}
	}

	/**
	 * Returns a vector of the first spatial derivatives in x, y & z 
	 * (nabla cC - see http://en.wikipedia.org/wiki/Del).
	 * Does this by first converting the ContinuousVector to a DiscreteVector and
	 * then estimating then gradient using the Mean Value Theorem 
	 * (http://en.wikipedia.org/wiki/Mean_value_theorem).
	 * @param cC
	 * @return
	 */
	public ContinuousVector getGradient(ContinuousVector cC) {
		DiscreteVector dV = new DiscreteVector(cC,_reso);
		return new ContinuousVector(diffX(dV), diffY(dV), diffZ(dV));
	}

	public ContinuousVector getGradient2D(ContinuousVector cC) {
		DiscreteVector dV = new DiscreteVector(cC,_reso);
		return new ContinuousVector(diffX(dV), diffY(dV), diffY(dV));
	}

	public ContinuousVector getGradient(DiscreteVector dV) {
		return new ContinuousVector(diffX(dV), diffY(dV), diffZ(dV));
	}

	/**
	 * Send the value on the padded grid at a given position
	 */
	public double getValueAt(DiscreteVector dV) {
		//sonia:chemostat

		if(Simulator.isChemostat){
			if (isValid(dV)) return grid[dV.i][dV.j][dV.k];
			else return Double.NaN;	
		}else{
			if (isValid(dV)) return grid[dV.i+1][dV.j+1][dV.k+1];
			else return Double.NaN;
		}
	}

	/**
	 * Send the value on the padded grid at a given position (the coordinates
	 * are NOT corrected)
	 */
	public double getValueAt(int i, int j, int k) {
		if (isValidorPadded(i, j, k)) return grid[i][j][k];
		else return Double.NaN;
	}

	/**
	 * Set a value on the padded grid (the coordinates are corrected)
	 * 
	 * @param value
	 * @param cc
	 */
	public void setValueAt(double value, ContinuousVector cC) {
		DiscreteVector dC = getDiscreteCoordinates(cC);
		setValueAt(value, dC);
	}

	/**
	 * Set a value on the padded grid (the coordinates are corrected)
	 * 
	 * @param value
	 * @param cc
	 */
	public void setValueAt(double value, DiscreteVector dC) {
		//sonia:chemostat
		if(Simulator.isChemostat){
			grid[dC.i][dC.j][dC.k] = value;
		}else{
			grid[dC.i+1][dC.j+1][dC.k+1] = value;
		}
	}

	/**
	 * Set a value on the padded grid (the coordinates are NOT corrected)
	 * 
	 * @param value
	 * @param i,j,k
	 */
	public void setValueAt(double value, int i, int j, int k) {
		grid[i][j][k] = value;
	}

	/**
	 * Set a value on the padded grid (the coordinates are corrected)
	 * 
	 * @param value
	 * @param cC
	 */
	public void addValueAt(double value, ContinuousVector cC) {
		addValueAt(value, getDiscreteCoordinates(cC));
	}

	/**
	 * Set a value on the padded grid (the coordinates are corrected)
	 * 
	 * @param value
	 * @param dC
	 */
	public void addValueAt(double value, DiscreteVector dC) {
		//sonia:chemostat
		if(Simulator.isChemostat){
			grid[dC.i][dC.j][dC.k] += value;
		}else{
			grid[dC.i+1][dC.j+1][dC.k+1] += value;
		}
	}

	/**
	 * Add a value on the padded grid (the coordinates are NOT corrected)
	 * 
	 * @param value
	 * @param i,j,k
	 */
	public void addValueAt(double value, int i, int j, int k) {
		grid[i][j][k] += value;
	}

	/**
	 * Add value to all grid locations (including the padding)
	 * 
	 * @param value
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
	 * Set to zero negative values in the grid
	 * 
	 * @param i
	 * @param j
	 * @param k
	 */
	public void truncateValueAt(int i, int j, int k) {
		grid[i][j][k] = (grid[i][j][k]<0 ? 0 : grid[i][j][k]);
	}

	/**
	 * Set all meshes of a grid with the same value (including the padding)
	 * 
	 * @param value
	 */
	public void setAllValueAt(double value) {

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
	 * @return number of real voxels along X (ignore the padding)
	 */
	public int getGridSizeI() {
		return _nI;
	}

	/**
	 * @return number of real voxels along Y (ignore the padding)
	 */
	public int getGridSizeJ() {
		return _nJ;
	}

	/**
	 * @return number of real voxels along Z (ignore the padding)
	 */
	public int getGridSizeK() {
		return _nK;
	}

	/**
	 * @param axeCode : 1-X, 2-Y,3-Z
	 * @return : the number of voxels along a direction including padding bands
	 */
	public int getGridTotalSize(int axeCode) {
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
	 * @param axeCode : 1-X, 2-Y, 3-Z
	 * @return : the length (in distance unit) along a direction ignoring
	 * padding bands
	 */
	public double getGridLength(int axeCode) {
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

	public double getVoxelVolume() {
		return ExtraMath.cube(_reso);
	}

	/**
	 * @return the whole grid including padding band
	 */
	public double[][][] getGrid() {
		return grid;
	}

	public double[][][] getCloneGrid(){
		return grid.clone();
	}

	public double getResolution() {
		return _reso;
	}

	public boolean is3D() {
		//sonia:chemostat
		if(Simulator.isChemostat){
			return false;
		}else{
			return _is3D;
		}
	}

	public void setGrid(double[][][] u) {
		utils.MatrixOperations.copyValuesTo(grid,u);
	}

	public void setGrid(double[][][] u, double v) {
		utils.MatrixOperations.muliplyBy(u, v);
		grid = u;
	}

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

			// Fill the mark-up
			if (_nK==1) {
				// We have a 2D grid
				for (int i = 0; i<_nI+2; i++) {
					for (int j = 0; j<_nJ+2; j++) {
						//bufferState.write(Arrays.toString(grid[i][j]));
						bufferState.write(Double.toString(grid[i][j][1]));
						bufferState.write(";\n");
					}
				}
			} else {
				// We have a 3D grid 
				for (int i = 0; i<_nI+2; i++) {
					for (int j = 0; j<_nJ+2; j++) {
						for (int k=0; k<_nK+2; k++) {
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
