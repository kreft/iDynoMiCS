/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * Implements connected volume filtration operation. Provides a base class for
 * geometry specific connected volume filtrators.
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */


/*
 * Created on 29-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */

package simulator.detachment;

import simulator.agent.LocatedGroup;

public class ConnectedVolume {

	protected int            _nTotal;
	protected int[]          gridDim = new int[3];
	protected boolean[]      _cvf;
	protected LocatedGroup[] _shoveGrid;
	protected int            _validateInThisIteration;

	/**
	 * Initialise the _cvf matrix
	 */
	public ConnectedVolume(int nI, int nJ, int nK) {
		gridDim[0] = nI;
		gridDim[1] = nJ;
		gridDim[2] = nK;
		_nTotal = (nI+2)*(nJ+2)*(nK+2);
		_cvf = new boolean[_nTotal];
	}

	/**
	 * Implements iterative implementation of connected volume filtering
	 * Iteration stop when none grid cell has been added
	 * @param matrixToFilter : vectorised array of space occupation matrix
	 * @return a vectorized array of a boolean matrix (true for connected to
	 * carrier)
	 */
	public boolean[] computeCvf(LocatedGroup[] matrixToFilter) {
		// assign
		_shoveGrid = matrixToFilter;

		// initiate the result matrix (true for carrier)
		initializeCvf();

		// Iteration
		_validateInThisIteration = 1;
		while (_validateInThisIteration>0) {
			_validateInThisIteration = 0;
			for (int index = 0; index<_nTotal; index++) {
				if (validateElement(index)) _validateInThisIteration++;
			}
		}
		return _cvf;
	}

	/**
	 * Initiate the cvf process. Research carrier on the padded grid
	 */
	protected void initializeCvf() {
		int status;
		for (int index = 0; index<_nTotal; index++) {
			status = _shoveGrid[index].status;
			// First test if the element is outside a cyclic boundary
			if(status==-1) continue;
			// Now test if it is a carrier
			_cvf[index] = (status==0);
		}
	}

	/**
	 * @return true if this grid element has been validated
	 */
	public boolean validateElement(int index) {
		boolean test;

		if (!_cvf[index]&_shoveGrid[index].status==1) {
			// check if one of your neighbors is attached to the carrier
			// (Note that we only check the neighbors in the cube-face directions
			// and not along any of the diagonals.)

			test = _cvf[_shoveGrid[index].nbhIndex[1][1][1]];

			// top/bottom neighbors (X direction)
			test |= _cvf[_shoveGrid[index].nbhIndex[2][1][1]];
			test |= _cvf[_shoveGrid[index].nbhIndex[0][1][1]];

			// left/right neighbors (Y direction)
			test |= _cvf[_shoveGrid[index].nbhIndex[1][0][1]];
			test |= _cvf[_shoveGrid[index].nbhIndex[1][2][1]];

			if (gridDim[2]>1) {
				// if it is a 3d simulation, do the front/back neighbors (Z direction)
				test |= _cvf[_shoveGrid[index].nbhIndex[1][1][0]];
				test |= _cvf[_shoveGrid[index].nbhIndex[1][1][2]];
			}

			if (test) {
				_cvf[index] = true;
				return true;
			}
		}
		return false;
	}

	/**
	 * Implements iterative implementation of connected volume filtering, but does
	 * so using an iterative process starting from the carrier.
	 * 
	 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu)
	 * 
	 * @param matrixToFilter : vectorised array of space occupation matrix
	 * @return a vectorized array of a boolean matrix (true for connected to
	 * carrier)
	 * 
	 * @deprecated This was more of a test than meant to be used for all cases
	 */
	public boolean[] computeCvfFromCarrier(LocatedGroup[] matrixToFilter) {
		// assign
		_shoveGrid = matrixToFilter;

		// go through and initialize the connections by starting at carrier indices
		int status;
		for (int index = 0; index < _nTotal; index++) {
			status = _shoveGrid[index].status;
			// First test if the element is outside a boundary
			if (status==-1) continue;

			// Now test if it is a carrier; if so, connect elements out from here
			if (status == 0) {
				_cvf[index] = true;
				attachElementNeighbors(index);
			}
		}

		return _cvf;
	}

	/**
	 * This is called from an attached element and looks at whether any
	 * neighbors should also be attached.
	 * This is first called by computeCvfFromCarrier().
	 * 
	 * @deprecated This was more of a test than meant to be used for all cases
	 */
	public void attachElementNeighbors(int index) {
		int testindex;

		// iterate over neighbors and attach any that have biomass AND aren't
		// already attached
		if (gridDim[2] > 1) {
			// 3D geometry
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					for (int k=0; k<3; k++) {
						testindex = _shoveGrid[index].nbhIndex[i][j][k];
						if (!_cvf[testindex] & _shoveGrid[testindex].status==1) {
							_cvf[testindex] = true;
							attachElementNeighbors(testindex);
						}
					}
		} else {
			// 2D geometry
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++) {
					testindex = _shoveGrid[index].nbhIndex[i][j][1];
					if (!_cvf[testindex] & _shoveGrid[testindex].status==1) {
						_cvf[testindex] = true;
						attachElementNeighbors(testindex);
					}
				}
		}
	}
}