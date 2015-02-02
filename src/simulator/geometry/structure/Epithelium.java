package simulator.geometry.structure;

import simulator.agent.zoo.EpithelialCell;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.geometry.boundaryConditions.BoundaryEpithelium;
import simulator.geometry.pointProcess.Site;
import simulator.geometry.pointProcess.Voronoi;
import utils.LogFile;

public abstract class Epithelium
{
	private BoundaryEpithelium[] _internalBCs;
	
	private AllBC[] _externalBCs;
	
	private EpithelialCell[] _cells;
	
	private Voronoi voronoi;
	
	public void setBoundaries(BoundaryEpithelium[] internal, AllBC[] external)
	{
		_internalBCs = internal;
		_externalBCs = external;
	}
	
	public void setCells(EpithelialCell[] cells)
	{
		_cells = cells;
	}
	
	public void solveVoronoi()
	{
		assignCells();
		
		for ( BoundaryEpithelium bEpi : _internalBCs )
		{
			voronoi = new Voronoi(bEpi.getShape(), bEpi.sites);
		}
	}
	
	private void assignCells()
	{
		Site cellSite;
		for ( BoundaryEpithelium bEpi : _internalBCs )
			bEpi.sites.clear();
		
		cellLoop: for ( EpithelialCell cell : _cells )
		{
			cellSite = new Site(cell.getLocation());
			for ( BoundaryEpithelium bEpi : _internalBCs )
				if ( bEpi.getDistance(cell.getLocation()) < 1E-6 )
				{
					bEpi.sites.add(cellSite);
					continue cellLoop;
				}
			LogFile.writeLogAlways("Epithelium could not place cell at "+
											cell.getLocation().toString());
		}
	}
}
