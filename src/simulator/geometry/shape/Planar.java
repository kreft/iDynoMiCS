
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package simulator.geometry.shape;

import simulator.geometry.*;
import simulator.SpatialGrid;
import utils.LogFile;
import utils.XMLParser;

import java.io.Serializable;
import java.util.*;

public class Planar implements IsShape, Serializable {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	/* Definition of the plan _______________________________________ */

	// Computation domain where this shape is defined

	// A point on the plan and a vector normal to the plan and going outside the
	// domain
	private DiscreteVector          _pointDCIn, _vectorDCOut;
	private ContinuousVector        _pointIn, _vectorOut;

	// Two orthogonal vectors colinear to the plan
	private DiscreteVector          u, v;
	// Range of discrete coordinates met on this shape
	private int                     uMax, vMax;

	/* Temporary variables __________________________________________ */
	private static ContinuousVector tempVar = new ContinuousVector();
	private static int              indexU, indexV;
	private static DiscreteVector   move    = new DiscreteVector();
	private static DiscreteVector   origin  = new DiscreteVector();

	/**
     * Tune a shape
     * @param shapeRoot
     * @param aDomain
     */
	public void readShape(XMLParser shapeRoot, Domain aDomain) {

		// Build the variables describing the plan
		_pointDCIn = shapeRoot.getParamXYZ("pointIn");
		_vectorDCOut = shapeRoot.getParamXYZ("vectorOut");

		// Translate them into continuous coordinates
		double res = aDomain.getGrid().getResolution();
		_pointIn = new ContinuousVector();
		_pointIn.x = (_pointDCIn.i+(1-_vectorDCOut.i)/2)*res;
		_pointIn.y = (_pointDCIn.j+(1-_vectorDCOut.j)/2)*res;
		_pointIn.z = (_pointDCIn.k+(1-_vectorDCOut.k)/2)*res;

		_vectorOut = new ContinuousVector();
		_vectorOut.x = _vectorDCOut.i;
		_vectorOut.y = _vectorDCOut.j;
		_vectorOut.z = _vectorDCOut.k;

		// Find two orthogonal vectors colinear to the plan
		u = new DiscreteVector();
		v = new DiscreteVector();
		_vectorDCOut.orthoVector(u, v);
	}

	/**
     * Test if the given coordinates are outside the boundary
     */
	public Boolean isOutside(ContinuousVector cc) {
		tempVar.x = -_pointIn.x+cc.x;
		tempVar.y = -_pointIn.y+cc.y;
		tempVar.z = -_pointIn.z+cc.z;

		return (_vectorOut.cosAngle(tempVar)>0);
	}

	/**
     * Computes orthogonal distance and if this distance is lower than the
     * resolution and if the point is outside, then it is declared on boundary
     */
	public Boolean isOnBoundary(ContinuousVector cC, double res) {
		return (isOutside(cC)&&(cC.distance(getOrthoProj(cC))<=res));
	}

	/**
     * @return : coordinates of the intersection between a line (described by a
     * point and a vector ) and the plane ; return null if none intersection
     * exists
     */
	public ContinuousVector intersection(ContinuousVector position, ContinuousVector vector) {

		// Determine the constant term for the equation of the plane
		double d = -_vectorOut.prodScalar(_pointIn);
		if (_vectorOut.prodScalar(vector)==0) {
			// the line will never cross this plane
			return null;
		}

		double k = (-d-_vectorOut.prodScalar(position))/_vectorOut.prodScalar(vector);

		ContinuousVector out = new ContinuousVector();
		out.x = position.x+k*vector.x;
		out.y = position.y+k*vector.y;
		out.z = position.z+k*vector.z;
		return out;
	}

	/**
     * Send a normal vector pointing toward the inside
     */
	public ContinuousVector getNormalInside(ContinuousVector cc) {
		return new ContinuousVector(-_vectorOut.x, -_vectorOut.y, -_vectorOut.z);
	}

	public ContinuousVector getNormalOutside(ContinuousVector cc) {
		return new ContinuousVector(_vectorOut);
	}

	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut) {
		double a, b, c, d, k;
		a = _vectorOut.x;
		b = _vectorOut.y;
		c = _vectorOut.z;
		d = -(_pointIn.x*a+_pointIn.y*b+_pointIn.z*c);

		// this next line wasn't calculating the projection coefficient correctly
		// and would case trouble when the plane was away from the substratum
		//k = (d+_vectorOut.prodScalar(ccIn))/(a+b+c);
		// this does it right
		k = -(d+_vectorOut.prodScalar(ccIn))/_vectorOut.prodScalar(_vectorOut);
		ccOut.x = ccIn.x+k*a;
		ccOut.y = ccIn.y+k*b;
		ccOut.z = ccIn.z+k*c;
	}

	public ContinuousVector getOrthoProj(ContinuousVector ccIn) {
		ContinuousVector ccOut = new ContinuousVector();
		orthoProj(ccIn, ccOut);
		return ccOut;
	}

	public double getDistance(IsShape aShape) {		
		ContinuousVector ccOut = aShape.intersection(_pointIn, _vectorOut);
		return _pointIn.distance(ccOut);
	}
	
	

	public double getDistance(ContinuousVector cc){
		ContinuousVector ccOut = intersection(cc, _vectorOut);
		return ccOut.distance(cc);
	}

	/**
     * Initialisation to go along the boundary
     */
	public void readyToFollowBoundary(SpatialGrid aSG) {
		double res = aSG.getResolution();

		origin.i = ((int)Math.floor(_pointIn.x/res))-u.i;
		origin.j = ((int)Math.floor(_pointIn.y/res))-u.j;
		origin.k = ((int)Math.floor(_pointIn.z/res))-u.k;
		
		if(_vectorDCOut.i>0) origin.i+=-1;
		if(_vectorDCOut.j>0) origin.j+=-1;
		if(_vectorDCOut.k>0) origin.k+=-1;

		indexU = 0;
		indexV = 0;
		
		if (u.i!=0) uMax = aSG.getGridSizeI()/u.i;
		else uMax = 0;
		if (u.j!=0) uMax = Math.max(uMax, aSG.getGridSizeJ()/u.j);
		if (u.k!=0) uMax = Math.max(uMax, aSG.getGridSizeK()/u.k);

		if (v.i!=0) vMax = aSG.getGridSizeI()/v.i;
		else vMax = 0;
		if (v.j!=0) vMax = Math.max(vMax, aSG.getGridSizeJ()/v.j);
		if (v.k!=0) vMax = Math.max(vMax, aSG.getGridSizeK()/v.k);
	}

	/**
     * 
     */
	public boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut, SpatialGrid aSG) {
	    // Find the next valid point
	    boolean vectorValid = false;
	    do {
		stepBoundary();
		dcIn.sendSum(origin, move);
		vectorValid = aSG.isValid(dcIn);
		// If a valid point has been found, compute its closest neighbour outside
		if (vectorValid) dcOut.sendSum(dcIn, _vectorDCOut);
	    } while ( !(vectorValid) && indexV<vMax );

	    return vectorValid;
	}

	/**
     * 
     */
	public void stepBoundary() {
		if (indexU<uMax) indexU++;
		else {
			indexU = 0;
			indexV++;
		}
		move.i = indexU*u.i+indexV*v.i;
		move.j = indexU*u.j+indexV*v.j;
		move.k = indexU*u.k+indexV*v.k;
	}

	@Deprecated
	public LinkedList<DiscreteVector> sendOnBorder(SpatialGrid aSoluteGrid) {

		double res = aSoluteGrid.getResolution();

		double x, y, z;
		double a, b, c, d;
		a = _vectorOut.x;
		b = _vectorOut.y;
		c = _vectorOut.z;
		d = -(_pointIn.x*a+_pointIn.y*b+_pointIn.z*c);

		// infinetisimal mouvement toward the outside
		double epsilon = 1e-6;

		LinkedList<DiscreteVector> out = new LinkedList<DiscreteVector>();
		ContinuousVector cC;
		if (b!=0) {

			for (int k = -1; k<=aSoluteGrid.getGridSizeK(); k++) {
				for (int i = -1; i<=aSoluteGrid.getGridSizeI(); i++) {
					y = (a*(i*res)+c*(k*res)+d)/b;
					cC = new ContinuousVector(i*res, y, k*res);
					cC.x += epsilon/2*_vectorOut.x;
					cC.y += epsilon/2*_vectorOut.y;
					cC.z += epsilon/2*_vectorOut.z;
					if (aSoluteGrid.isValidOrPadded(cC)&&isOnShape(cC)) {
						out.add(aSoluteGrid.getDiscreteCoordinates(cC));
					}
				}
			}
			return out;
		}

		if (a!=0) {
			for (int j = -1; j<=aSoluteGrid.getGridSizeJ(); j++) {
				for (int k = -1; k<=aSoluteGrid.getGridSizeK(); k++) {
					x = (b*(j*res)+c*(k*res)+d)/a;
					cC = new ContinuousVector(x, j*res, k*res);
					cC.x += epsilon/2*_vectorOut.x;
					cC.y += epsilon/2*_vectorOut.y;
					cC.z += epsilon/2*_vectorOut.z;
					if (aSoluteGrid.isValidOrPadded(cC)) {
						out.add(aSoluteGrid.getDiscreteCoordinates(cC));
					}
				}
			}
			return out;
		}

		if (c!=0) {
			for (int j = -1; j<=aSoluteGrid.getGridSizeJ(); j++) {
				for (int i = -1; i<=aSoluteGrid.getGridSizeI(); i++) {
					z = (a*(i*res)+b*(j*res)+d)/c;
					cC = new ContinuousVector(i*res, j*res, z);
					cC.x += epsilon/2*_vectorOut.x;
					cC.y += epsilon/2*_vectorOut.y;
					cC.z += epsilon/2*_vectorOut.z;
					if (aSoluteGrid.isValidOrPadded(cC)) {
						out.add(aSoluteGrid.getDiscreteCoordinates(cC));
					}
				}
			}
			return out;
		}
		return out;
	}

	@Deprecated
	public LinkedList<ContinuousVector> sendCCOnBorder(SpatialGrid aSoluteGrid) {

		double res = aSoluteGrid.getResolution();

		LinkedList<ContinuousVector> out = new LinkedList<ContinuousVector>();
		ContinuousVector cC;
		int iMin = -(int) Math.floor(_pointIn.x/res)-1;
		int iMax = (int) Math.floor(aSoluteGrid.getGridSizeI()-_pointIn.x/res)+1;
		int jMin = -(int) Math.floor(_pointIn.y/res)-1;
		int jMax = (int) Math.floor(aSoluteGrid.getGridSizeJ()-_pointIn.y/res)+1;
		int kMin = -(int) Math.floor(_pointIn.z/res)-1;
		int kMax = (int) Math.floor(aSoluteGrid.getGridSizeK()-_pointIn.z/res)+1;

		ContinuousVector nV1 = new ContinuousVector();
		nV1.x = Math.abs(_vectorOut.z);
		nV1.y = Math.abs(_vectorOut.x);
		nV1.z = Math.abs(_vectorOut.y);

		ContinuousVector nV2 = new ContinuousVector();
		nV2.x = nV1.z;
		nV2.y = nV1.x;
		nV2.z = nV1.y;

		for (int i = iMin; i<=iMax; i++) {
			for (int j = jMin; j<=jMax; j++) {
				for (int k = kMin; k<=kMax; k++) {
					cC = new ContinuousVector(_pointIn);
					cC.x += i*(nV1.x+nV2.x)*res;
					cC.y += j*(nV1.y+nV2.y)*res;
					cC.z += k*(nV1.z+nV2.z)*res;
					if (aSoluteGrid.isValidOrPadded(cC)) {
						out.add(cC);
					}
				}
			}

		}

		return out;
	}

	/**
     * test that the received point is coplanar with the definition point
     * @param cC
     * @return
     */
	public boolean isOnShape(ContinuousVector cC) {
		boolean out = ((cC.x-_pointIn.x)/_vectorOut.x==(cC.y-_pointIn.y)/_vectorOut.y);
		out &= ((cC.x-_pointIn.x)/_vectorOut.x==(cC.z-_pointIn.z)/_vectorOut.z);
		return out;
	}

	public ContinuousVector getPointIn() {
		return _pointIn;
	}

	public DiscreteVector getNormalDC() {
		return _vectorDCOut;
	}
}
