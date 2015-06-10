package odeSolver;

import idyno.SimTimer;

import java.util.ArrayList;

import simulator.reaction.Reaction;
import utils.LogFile;
import Jama.Matrix;

public class ChemostatSolver extends ODEsolver
{
	protected Double _dilution;
	
	protected ArrayList<Reaction> _reactions;
	
	protected int _nReac;
	
	protected Matrix _inflow;
	
	protected Double[] _catalysts;
	
	public void init(int nVar, Double dilution, ArrayList<Reaction> reactions,
													Double hmax, Double rtol)
	{
		super.init(nVar, false, hmax, rtol);
		this._dilution = dilution;
		this._reactions = reactions;
		this._nReac = this._reactions.size();
	}
	
	public void updateInflow(Matrix inflow)
	{
		this._inflow = inflow;
	}
	
	public Matrix solve(Matrix solutes, Double[] catalysts)
	{
		this._catalysts = catalysts;
		return this.solve(solutes, SimTimer.getCurrentTimeStep());
	}
	
	/**
	 * 
	 * @param y
	 * @param target 
	 */
	protected Matrix calc1stDeriv(Matrix y, Matrix target)
	{
		target = this._inflow.copy();
		target.minusEquals(y);
		target.timesEquals(this._dilution);
		for ( int iReac = 0; iReac < this._nReac; iReac++ )
		{
			target.plusEquals(
				this._reactions.get(iReac).calcdMUdT(y, this._catalysts[iReac]));
		}
		return target;
	}
	
	/**
	 * Update the second derivative of Y, i.e. the rate of change of F with
	 * respect to time (dFdT).
	 * 
	 * @param y 
	 * @param tdel
	 * @param target  
	 */
	protected Matrix calc2ndDeriv(Matrix y, Double tdel, Matrix target)
	{
		return this.numerical2ndDeriv(y, tdel, target);
	}
	
	/**
	 * Update the Jacobian matrix, i.e. the rate of change of F with respect to
	 * each of the variables in Y (dFdY).
	 * 
	 * @param y 
	 * @param target 
	 */
	protected Matrix calcJacobian(Matrix y, Matrix target)
	{
		target = Matrix.identity(this._nVar, this._nVar).times(-this._dilution);
		Matrix dMUdS;
		for ( int iReac = 0; iReac < this._nReac; iReac++ )
		{
			dMUdS = this._reactions.get(iReac).calcdMUdS(y, this._catalysts[iReac]);
			target.plusEquals(dMUdS);
		}
		return target;
	}
}