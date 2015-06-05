package odeSolver;

import Jama.Matrix;
import simulator.agent.zoo.GeneRegBac;

/**
 * 
 * 
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class GeneRegSolver extends ODEsolver
{
	/**
	 * The cell to be referred to when calculating the necessary derivatives.
	 * 
	 * @see 
	 */
	protected GeneRegBac _refCell;
	
	
	public GeneRegSolver()
	{
		super();
	}
	
	/**
	 * 
	 */
	public void init(int nVar)
	{
		super.init(nVar, false);
		//super.init(nVar, true);
	}
	
	public void setReferenceAgent(GeneRegBac aCell)
	{
		_refCell = aCell;
	}
	
	/**
	 * 
	 */
	public Matrix calc1stDeriv(Matrix y)
	{
		return _refCell.calc1stDeriv(y);
	}
	
	
	/**
	 * Update the second derivative of Y, i.e. the rate of change of F with
	 * respect to time (dFdT).
	 * 
	 * It is implicitly assumed that solute concentrations are constant
	 * throughout the gene regulation solution, i.e. that gene regulation
	 * happens on a much faster timescale than the solute dynamics.
	 * 
	 * @param y 
	 * @param tdel 
	 */
	public Matrix calc2ndDeriv(Matrix y, Double tdel)
	{
		Matrix dFdT = new Matrix(nVar, 1, 0.0);
		
		Matrix dYdT = calc1stDeriv(y);
		
		// yNext = levels + (tdel * dYdT)
		Matrix yNext = y.plus(dYdT.times(tdel));
		
		// dFdT = ( dYdT(ynext) - dYdT(y) )/tdel
		dFdT = calc1stDeriv(yNext);
		dFdT.minusEquals(dYdT);
		dFdT.timesEquals(1/tdel);
		
		return dFdT;
	}
	
	/**
	 * Update the Jacobian matrix, i.e. the rate of change of F with respect to
	 * each of the variables in Y (dFdY).
	 * 
	 * @param levels Matrix (column) of the gene expression levels.
	 * @return 
	 */
	public Matrix calcJacobian(Matrix levels)
	{
		return _refCell.calcJacobian(levels);
	}
}
