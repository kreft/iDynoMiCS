/**
 * \package simulator.geometry.shape
 * \brief Package of utilities that assist with managing the boundary conditions in iDynoMiCS
 * 
 * Package of utilities that assist with managing the boundary conditions in iDynoMiCS. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.geometry.shape;

import simulator.geometry.*;
import simulator.SpatialGrid;
import utils.XMLParser;
import java.io.Serializable;
import java.util.*;
import utils.LogFile;

/**
 * \brief Create a planar shaped boundary
 * 
 * Each boundaryCondition also includes a shape mark-up to define the shape of the boundary. In this release only the 'Planar' class is 
 * available, and it requires specification of a point on the boundary and a vector pointing outside the domain. These shape parameters 
 * must be given in index coordinates 
 * 
 */
public class Planar implements IsShape, Serializable 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;
	/* Definition of the plan _______________________________________ */

	// Computation domain where this shape is defined

	/**
	 * A point on the plane 
	 */ 
	private DiscreteVector          _pointDCIn;
	
	/**
	 * A vector normal to the plan and going outside the domain
	 */
	private DiscreteVector			_vectorDCOut;
	
	/**
	 * A point on the plane (ContinuousVector)
	 */ 
	private ContinuousVector        _pointIn;
	
	/**
	 * A vector normal to the plan and going outside the domain (ContinuousVector)
	 */
	private ContinuousVector		 _vectorOut;

	/**
	 *  Orthogonal vectors colinear to the plan
	 */
	private DiscreteVector          u;
	
	/**
	 *  Orthogonal vectors colinear to the plan
	 */
	private DiscreteVector			v;
	
	/**
	 * Range of discrete coordinates met on this shape
	 */
	private int                     uMax;
	
	/**
	 * Range of discrete coordinates met on this shape
	 */
	private int						vMax;

	/* Temporary variables __________________________________________ */
	
	/**
	 * Temporary continuous vector to use in calculations
	 */
	private static ContinuousVector tempVar = new ContinuousVector();
	
	/**
	 * Index used to check whether a point is within the shape
	 */
	private static int              indexU;
	
	/**
	 * Index used to check whether a point is within the shape
	 */
	private static int				indexV;
	
	/**
	 * Stores the move while being calculated
	 */
	private static DiscreteVector   move    = new DiscreteVector();
	
	/**
	 * Temporary store of the origin of a point prior to move
	 */
	private static DiscreteVector   origin  = new DiscreteVector();

	/**
	 * \brief Reads the coordinates that specify a boundary from the protocol file, creating a shape
	 * 
	 * Reads the coordinates that specify a boundary from the protocol file, creating a shape
	 * 
	 * @param shapeRoot	XML elements from the protocol file that contain coordinates specifying the edge of a boundary
	 * @param aDomain	The computation domain that this boundary is associated with
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain) 
	{

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
	 * \brief Test if a given set of coordinates is outside this shape
	 * 
	 * Test if a given set of coordinates is outside this shape
	 * 
	 * @param cc	ContinuousVector containing the coordinates of a point to test
	 * @return	Boolean noting whether this coordinate is inside or outside this shape
	 */
	public Boolean isOutside(ContinuousVector cc) {
		tempVar.x = -_pointIn.x+cc.x;
		tempVar.y = -_pointIn.y+cc.y;
		tempVar.z = -_pointIn.z+cc.z;

		return (_vectorOut.cosAngle(tempVar)>0);
	}

	/**
     * Computes orthogonal distance and if this distance is lower than the resolution and if the point is outside, 
     * then the point tested is declared to be on the boundary of the domain
     * 
	 * @param cC	ContinuousVector containing the coordinates of a point to test
	 * @param res	Resolution of the domain that this shape is associated with
	 * @return	Boolean noting whether this coordinate is on the boundary of the domain
	 */
	public Boolean isOnBoundary(ContinuousVector cC, double res) {
		return (isOutside(cC)&&(cC.distance(getOrthoProj(cC))<=res));
	}

	/**
     * \brief Calculates the coordinates of the interaction between a line (point a vector) and the plane
     * 
     * Calculates the coordinates of the interaction between a line (point a vector) and the plane. Returns null if none exists
     * 
     * @param position	Position used to calculate the line
     * @param vector	Vector of coordinate positions used to calculate the line
     * @return : coordinates of the intersection between a line and the plane
     */
	public ContinuousVector intersection(ContinuousVector position, ContinuousVector vector) {

		// Determine the constant term for the equation of the plane
		double d = _vectorOut.prodScalar(_pointIn);
		if (_vectorOut.prodScalar(vector)==0) {
			// the line will never cross this plane
			return null;
		}

		double k = (d-_vectorOut.prodScalar(position))/_vectorOut.prodScalar(vector);

		ContinuousVector out = new ContinuousVector();
		out.x = position.x+k*vector.x;
		out.y = position.y+k*vector.y;
		out.z = position.z+k*vector.z;
		return out;
	}

	/**
	 * \brief Takes a vector and returns that vector pointing towards the inside of the shape
	 * 
	 * Takes a vector and returns that vector pointing towards the inside of the shape
	 * 
	 * @param cc	Vector outside the shape
	 * @return ContinuousVector that is pointing towards the inside of the shape
	 * 
	 */
	public ContinuousVector getNormalInside(ContinuousVector cc) {
		return new ContinuousVector(-_vectorOut.x, -_vectorOut.y, -_vectorOut.z);
	}

	/**
	 * \brief Correct coordinates of a point that has gone outside this shape 
	 * 
	 * Correct coordinates of a point that has gone outside this shape
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @param ccOut	Corrected coordinates
	 * 
	 */
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

	/**
	 * \brief Correct coordinates of a point that has gone outside this shape, returning these coordinates
	 * 
	 * Correct coordinates of a point that has gone outside this shape, returning these coordinates
	 * 
	 * @param ccIn	Coordinates to be corrected
	 * @return Corrected coordinates
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn) {
		ContinuousVector ccOut = new ContinuousVector();
		orthoProj(ccIn, ccOut);
		return ccOut;
	}

	/**
	 * \brief Used in cyclic boundaries - gets the distance from the opposite side (aShape)
	 * 
	 * Used in cyclic boundaries - gets the distance from the opposite side (aShape)
	 * 
	 * @return Double stating distance to that shape
	 */
	public double getDistance(IsShape aShape) {		
		ContinuousVector ccOut = aShape.intersection(_pointIn, _vectorOut);
		return _pointIn.distance(ccOut);
	}
	
	
	/**
	 * \brief Used in cyclic boundaries - gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * Used in cyclic boundaries - gets the distance from a point on the other side (ContinuousVector)
	 * 
	 * @return Double stating distance to that shape
	 */
	public double getDistance(ContinuousVector cc){
		ContinuousVector ccOut = intersection(cc, _vectorOut);
		return ccOut.distance(cc);
	}

	/**
     * \brief Initialisation to create the features of and go along the boundary
     * 
     * Initialisation to create the features of and go along the boundary
     * 
     * @param aSG	The grid to which this boundary is a part
     */
	public void readyToFollowBoundary(SpatialGrid aSG) 
	{
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
     * \brief Find the next valid point
     * 
     *  Find the next valid point 
     *  
     *  @param dcIn	Discrete vector within the shape
     *  @param dcOut	Discrete vector outside the shape
     *  @param aSG	Spatial grid in which the boundary is associated
     *  @return Whether a valid point was found
     *  
     */
	public boolean followBoundary(DiscreteVector dcIn, DiscreteVector dcOut, SpatialGrid aSG) 
	{
		// Find the next valid point
		boolean vectorValid = false;
		do 
		{
			stepBoundary();
			dcIn.sendSum(origin, move);
			vectorValid = aSG.isValid(dcIn);
			// If a valid point has been found, compute its closest neighbour outside
			if (vectorValid) dcOut.sendSum(dcIn, _vectorDCOut);
		} while ( !(vectorValid) && indexV<vMax );
		
		return vectorValid; 
	}

	/**
     * \brief Process next location in the boundary
     * 
     * Process next location in the boundary
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

	/**
     * \brief Test that the received point is coplanar with the definition point
     * 
     * Test that the received point is coplanar with the definition point
     * 
     * @param cC	ContinuousVector of coordinates that should be tested
     * @return	Boolean stating whether the point is cooplanar with the definition point
     */
	public boolean isOnShape(ContinuousVector cC) {
		boolean out = ((cC.x-_pointIn.x)/_vectorOut.x==(cC.y-_pointIn.y)/_vectorOut.y);
		out &= ((cC.x-_pointIn.x)/_vectorOut.x==(cC.z-_pointIn.z)/_vectorOut.z);
		return out;
	}

	
	/**
	 * \brief Return vector normal to the plane
	 * 
	 * Return vector normal to the plane
	 * 
	 * @return	Discrete vector normal to the plane
	 */
	public DiscreteVector getNormalDC() {
		return _vectorDCOut;
	}
}
