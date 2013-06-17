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

import java.util.LinkedList;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.geometry.boundaryConditions.AllBC;

/**
 * \brief Interface describing mandatory methods for a computation domain object. The default class is "ComputationDomain" but other can be defined
 * 
 * Interface describing mandatory methods for a computation domain object. The default class is "ComputationDomain" but other can be defined
 * 
 * @author Andreas D�tsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public interface IsComputationDomain 
{

	/**
     * \brief Test if a given location is outside a boundary. Used to detect the crossed boundary when moving an agent
     * 
     * Test if a given location is outside a boundary. Used to detect the crossed boundary when moving an agent
     * 
     * @param newLoc	The location to test
     * @return	Boundary that the point has crossed (if applicable - null if no boundary crossed)
     */
	public AllBC testCrossedBoundary(ContinuousVector newLoc);

	/**
     * \brief Returns the domain grid
     * 
     * Returns the domain grid
     * 
     * @return SpatialGrid within this domain
     */
	public SpatialGrid getGrid();
	
	/**
	 * \brief Return the name of this domain
	 * 
	 * Return the name of this domain
	 * 
	 * @return Name of this domain
	 */
	public String getName();

	/**
     * \brief Return all boundaries associated with this domain
     * 
     * Return all boundaries associated with this domain
     * 
     * @return LinkedLIst of all boundaries in the container boundaryList
     */
	public LinkedList<AllBC> getAllBoundaries();
	
	
	/**
	 * \brief Return the boundary layer grid associated with this domain
	 * 
	 * Return the boundary layer associated with this domain
	 * 
	 * @return	SoluteGrid containing boundary between bulk and biofilm
	 */
	public SoluteGrid getBoundaryLayer();
	
	/**
	 * \brief Returns a list of discrete vectors that specify the limit of the boundary layer
	 * 
	 * Returns a list of discrete vectors that specify the limit of the boundary layer
	 * 
	 * @return LinkedList of DiscreteVectors with the limit of the boundary layer
	 */
	public LinkedList<DiscreteVector> getBorder();
	
	/**
	 * \brief Return the diffusivity grid associated with this domain
	 * 
	 * Return the diffusivity grid associated with this domain
	 * 
	 * @return	SoluteGrid containing diffusivity grid statistics
	 */
	public SoluteGrid getDiffusivity();
	
	/**
	 * \brief Return the biomass grid associated with this domain
	 * 
	 * Return the biomass grid associated with this domain
	 * 
	 * @return	SoluteGrid containing biomass throughout this domain
	 */
	public SoluteGrid getBiomass();
	
	/**
	 * \brief Refresh relative diffusivity and boundary layer grids to ensure biomass updated this step is included
	 * 
	 * Method to refresh relative diffusivity and boundary layer grids to ensure biomass updated this step is included. 
	 * Used in the creation of output files
	 */
	public void refreshBioFilmGrids();

}
