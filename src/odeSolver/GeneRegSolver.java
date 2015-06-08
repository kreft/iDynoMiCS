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
	}
	
	public void setReferenceAgent(GeneRegBac aCell)
	{
		_refCell = aCell;
	}
	
	/**
	 * 
	 * @param y
	 * @param target 
	 */
	@Override
	protected void calc1stDeriv(Matrix y, Matrix target)
	{
		target = _refCell.calc1stDeriv(y);
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
	 * @param target 
	 */
	protected void calc2ndDeriv(Matrix y, Double tdel, Matrix target)
	{
		this.numerical2ndDeriv(y, tdel, target);
	}
	
	/**
	 * Update the Jacobian matrix, i.e. the rate of change of F with respect to
	 * each of the variables in Y (dFdY).
	 * 
	 * @param levels Matrix (column) of the gene expression levels.
	 * @param target 
	 */
	protected void calcJacobian(Matrix levels, Matrix target)
	{
		target = _refCell.calcJacobian(levels);
	}
}
