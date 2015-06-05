package odeSolver;

import Jama.Matrix;

public class ChemostatSolver extends ODEsolver
{
	
	
	
	public Matrix calc1stDeriv(Matrix y)
	{
		Matrix out = new Matrix(nVar, 1, 0.0);
		// TODO
		return out;
	}
	
	/**
	 * Update the second derivative of Y, i.e. the rate of change of F with
	 * respect to time (dFdT).
	 * 
	 * @param y 
	 * @param tdel 
	 */
	public Matrix calc2ndDeriv(Matrix y, Double tdel)
	{
		Matrix out = new Matrix(nVar, 1, 0.0);
		// TODO
		return out;
	}
	
	/**
	 * Update the Jacobian matrix, i.e. the rate of change of F with respect to
	 * each of the variables in Y (dFdY).
	 * 
	 * @param y 
	 */
	public Matrix calcJacobian(Matrix y)
	{
		Matrix out = new Matrix(nVar, nVar, 0.0);
		// TODO
		return out;
	}
}