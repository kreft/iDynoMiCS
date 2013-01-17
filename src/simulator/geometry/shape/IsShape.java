
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package simulator.geometry.shape;

import java.util.LinkedList;

import simulator.geometry.Domain;
import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.SpatialGrid;
import utils.XMLParser;

public interface IsShape {

	public void readShape(XMLParser shapeRoot,Domain aDomain);
	public Boolean isOutside(ContinuousVector cc);
	public Boolean isOnBoundary(ContinuousVector cC, double res);	
	
	@Deprecated
	public LinkedList<DiscreteVector> sendOnBorder(SpatialGrid aSoluteGrid);
	@Deprecated
	public LinkedList<ContinuousVector> sendCCOnBorder(SpatialGrid aSoluteGrid);
	
	public ContinuousVector intersection(ContinuousVector position,ContinuousVector vector);
	public void orthoProj(ContinuousVector ccIn,ContinuousVector ccOut);
	public ContinuousVector getOrthoProj(ContinuousVector ccIn);
	
	public double getDistance(IsShape aShape);
	public double getDistance(ContinuousVector cc);
	public ContinuousVector getNormalInside(ContinuousVector cc);
	public ContinuousVector getNormalOutside(ContinuousVector cc);
	
	public void readyToFollowBoundary(SpatialGrid aSG);
	public boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut, SpatialGrid aSG);
	
	public ContinuousVector getPointIn();
	public DiscreteVector getNormalDC();

}
