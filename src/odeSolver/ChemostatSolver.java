package odeSolver;

import simulator.reaction.Reaction;
import Jama.Matrix;

public class ChemostatSolver extends ODEsolver
{
	protected Double _dilution;
	
	protected Reaction[] _reactions;
	
	protected int _nReac;
	
	protected Matrix _inflow;
	
	protected Double[] _catalysts;
	
	public void init(int nVar, Double dilution, Reaction[] reactions)
	{
		super.init(nVar, false);
		this._dilution = dilution;
		this._reactions = reactions;
		this._nReac = this._reactions.length;
	}
	
	public void solve(Matrix solutes, Double[] catalysts)
	{
		this._catalysts = catalysts;
		
		
	}
	
	/**
	 * 
	 * @param y
	 * @param target 
	 */
	protected void calc1stDeriv(Matrix y, Matrix target)
	{
		target = this._inflow.copy();
		target.minusEquals(y);
		target.timesEquals(this._dilution);
		for ( int iReac = 0; iReac < this._nReac; iReac++ )
		{
			target.plusEquals(
				this._reactions[iReac].calcdMUdT(y, this._catalysts[iReac]));
		}
	}
	
	/**
	 * Update the second derivative of Y, i.e. the rate of change of F with
	 * respect to time (dFdT).
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
	 * @param y 
	 * @param target 
	 */
	protected void calcJacobian(Matrix y, Matrix target)
	{
		target = Matrix.identity(this._nVar, this._nVar).times(-this._dilution);
		for ( int iReac = 0; iReac < this._nReac; iReac++ )
		{
			target.plusEquals(
				this._reactions[iReac].calcdMUdS(y, this._catalysts[iReac]));
		}
	}
}