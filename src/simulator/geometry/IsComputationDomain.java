/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *______________________________________________________
 * IsComputationDomain : interface describing mandatory methods for a 
 * computation domain object. The default class is "ComputationDomain" but other 
 * can be defined
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * ____________________________________________________________________________
 */

package simulator.geometry;

import java.util.LinkedList;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.geometry.boundaryConditions.AllBC;

public interface IsComputationDomain {

	/**
     * Test if a given location is outside a boundary ; used to detect the
     * crossed boundary when moving an agent
     * @param newLoc
     * @param radius
     * @param aBoundary : reference used to store the crossed boundary
     * @return true if a boundary has been crossed
     */
	public AllBC testCrossedBoundary(ContinuousVector newLoc);

	/**
     * @return the domain grid
     */
	public SpatialGrid getGrid();
	public String getName();

	/**
     * @param i
     * @return the ith boundary in the container boundaryList
     */
	public LinkedList<AllBC> getAllBoundaries();
	public AllBC getClosestBoundaryWithBulk(ContinuousVector cc);
	
	public SoluteGrid getBoundaryLayer();
	public LinkedList<DiscreteVector> getBorder();
	public SoluteGrid getDiffusivity();
	public SoluteGrid getBiomass();
	
	public void refreshBioFilmGrids();

}
