/**
 * \package simulator.detachment
 * \brief Package of classes that capture detachment of agents from the
 * biomass.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.detachment;

import simulator.agent.LocatedGroup;

/**
 * \brief Implements connected volume filtration operation.
 * 
 * Provides a base class for geometry specific connected volume filters.
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer
 * Center (NY, USA)
 *
 */
public class ConnectedVolume
{
	/**
	 * Total number of elements in the grid.
	 */
	protected int _nTotal;
	
	/**
	 * Whether or not to check in the Z direction.
	 */
	protected Boolean is3D;
	
	/**
	 * Connected volume filtration matrix of size _nTotal.
	 */
	protected Boolean[] _cvf;
	
	/**
	 * Shoving grid to use to determine which agents to mark for sloughing.
	 */
	protected LocatedGroup[] _shoveGrid;
	
	/**
	 * Integer noting whether a grid element needs to be validated (1) in an
	 * iteration or not (0).
	 */
	protected int _validateInThisIteration;

	/**
	 * \brief Initialise the _cvf matrix.
	 */
	public ConnectedVolume(int nI, int nJ, int nK)
	{
		is3D = (nK > 1);
		_nTotal = (nI+2)*(nJ+2)*(nK+2);
		_cvf = new Boolean[_nTotal];
	}

	/**
	 * \brief Implements iterative implementation of connected volume
	 * filtering.
	 * 
	 * Iteration stop when none grid cell has been added.
	 * 
	 * @param matrixToFilter	Vectorized array of space occupation matrix.
	 * @return a vectorized array of a boolean matrix (true for connected to
	 * carrier).
	 */
	public Boolean[] computeCvf(LocatedGroup[] matrixToFilter)
	{
		/*
		 * Assign.
		 */
		_shoveGrid = matrixToFilter;
		/*
		 * Initiate the result matrix (true for carrier).
		 */
		initializeCvf();
		/*
		 * Iteration.
		 */
		_validateInThisIteration = 1;
		while ( _validateInThisIteration > 0 )
		{
			_validateInThisIteration = 0;
			for ( int index = 0; index < _nTotal; index++ )
				if ( validateElement(index) )
					_validateInThisIteration++;
		}
		return _cvf;
	}
	
	/**
	 * \brief Initiate the connected volume filtration process.
	 * 
	 * Research carrier on the padded grid.
	 */
	protected void initializeCvf() 
	{
		int status;
		for ( int index = 0; index < _nTotal; index++ ) 
		{
			_cvf[index] = false;
			status = _shoveGrid[index].status;
			/*
			 * First test if the element is outside a cyclic boundary
			 */
			if ( status == -1 )
				continue;
			/*
			 * Now test if it is a carrier
			 */
			_cvf[index] = (status == 0);
		}
	}

	/**
	 * \brief Determine if this grid element has been validated.
	 * 
	 * @return true if this grid element has been validated.
	 */
	public boolean validateElement(int index)
	{
		boolean test;
		if ( !_cvf[index] && _shoveGrid[index].status == 1 ) 
		{
			/*
			 * Check if one of your neighbors is attached to the carrier. Note
			 * that we only check the neighbors in the cube-face directions
			 * and not along any of the diagonals.
			 */
			test = _cvf[_shoveGrid[index].nbhIndex[1][1][1]];
			/*
			 * Top/bottom neighbors (X direction).
			 */
			test = test || _cvf[_shoveGrid[index].nbhIndex[2][1][1]];
			test = test || _cvf[_shoveGrid[index].nbhIndex[0][1][1]];
			/*
			 * Left/right neighbors (Y direction).
			 */
			test = test || _cvf[_shoveGrid[index].nbhIndex[1][0][1]];
			test = test || _cvf[_shoveGrid[index].nbhIndex[1][2][1]];
			/*
			 * If it is a 3d simulation, do the front/back neighbors
			 * (Z direction).
			 */
			if ( is3D )
			{
				test = test || _cvf[_shoveGrid[index].nbhIndex[1][1][0]];
				test = test || _cvf[_shoveGrid[index].nbhIndex[1][1][2]];
			}
			
			if ( test )
				return _cvf[index] = true;
		}
		return false;
	}
}