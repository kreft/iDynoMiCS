package odeSolver;


import utils.LogFile;
import Jama.Matrix;

public abstract class ODEsolver
{
	/**
	 * 
	 */
	protected Double _d = 1.0 / (2.0 + Math.sqrt(2.0));
	
	/**
	 * 
	 */
	protected Double _e32  = 6.0 + Math.sqrt(2.0);
	
	/**
	 * The order of this method is 3.
	 */
	protected Double _power = 1.0/3.0;
	
	/**
	 * Numerical accuracy for EPS (error per step) 
	 */
	protected Double _sqrtE = Math.sqrt(2.22e-16);
	
	/**
	 * Error Per Step: the smallest positive floating-point number such that
	 * 1.0 + EPS > 1.0 
	 */
	protected Double _EPS = 2.22e-16;
	
	/**
	 * Number of variables in the solver.
	 */
	protected int _nVar;
	
	/**
	 * 
	 */
	protected Boolean _allowNegatives;
	
	protected Matrix _ynext;
	
	protected Matrix _dYdT;
	
	protected Matrix _dFdT;
	
	protected Matrix _dFdY;
	
	protected Matrix _f1, _f2, _W, _invW, _k1, _k2, _k3, _kaux;
	
	protected Double _t, _error, _tnext, _tdel, _h, _test;
	
	protected Boolean lastStep, noFailed, usingHMin;
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public ODEsolver()
	{
		
	}
	
	public void init(int nVar, Boolean allowNegatives)
	{
		this._nVar = nVar;
		this._allowNegatives = allowNegatives;
		
		_ynext = new Matrix(nVar, 1, 0.0);
		_dYdT  = new Matrix(nVar, 1, 0.0);
		_dFdT  = new Matrix(nVar, 1, 0.0);
		_dFdY  = new Matrix(nVar, nVar, 0.0);
		_f1    = new Matrix(nVar, 1, 0.0);
		_f2    = new Matrix(nVar, 1, 0.0);
		_W     = new Matrix(nVar, nVar, 0.0);
		_invW  = new Matrix(nVar, nVar, 0.0);
		_k1    = new Matrix(nVar, 1, 0.0);
		_k2    = new Matrix(nVar, 1, 0.0);
		_k3    = new Matrix(nVar, 1, 0.0);
		_kaux  = new Matrix(nVar, 1, 0.0);
	}
	
	/*************************************************************************
	 * KEY METHOS
	 ************************************************************************/
	
	/**
	 * 
	 * @param y
	 * @param tfinal
	 * @param rtol
	 * @param hmax
	 * @return
	 */
	public Matrix solve(Matrix y, Double tfinal, Double rtol, Double hmax)
	{
		/*
		 * First check that y is the correct size.
		 */
		if ( y.getRowDimension() != _nVar || y.getColumnDimension() != 1 )
			throw new IllegalArgumentException("Wrong matrix dimensions");
		/*
		 * Reset the matrices we will need.
		 */
		
		/*
		 * Control statement in case the maximum timestep size, hmax, is too
		 * large.
		 */
		while ( hmax > tfinal )
		{
			hmax *= 0.5;
			rtol *= 0.5;
		}
		/*
		 * First try a step size of hmax.
		 */
		_t = 0.0;
		lastStep  = false;
		_h = hmax;
		while ( ! lastStep )
		{
			/*
			 * If the next step gets us close to the end, we may as well
			 * just finish.
			 */
			if ( 1.05 * _h >= tfinal - _t )
			{
				_h = tfinal - _t;
				lastStep = true;
			}
			/*
			 * Update dFdT based on the mini-timestep tdel. The Jacobian
			 * matrix, dFdY, doesn't need this mini-timestep.
			 */
			_tdel = _sqrtE * (_t+_h);
			calc2ndDeriv(y, _tdel, _dFdT);
			calcJacobian(y, _dFdY);
			/*
			 * Try out this value of h, keeping a note of whether it ever
			 * fails.
			 */
			noFailed = true;
			usingHMin = false;
			while ( true )
			{
				_tnext = ( lastStep ) ? tfinal : _t + _h;
				/*
				 * The Rosenbrock method.
				 */
				try
				{
					/*
					 * W = I - h * d * dFdY
					 */
					_W = _dFdY.times( - _h * _d );
					_W.plusEquals( Matrix.identity(_nVar, _nVar) );
					if (_W.cond() > 10)
					{ 
						LogFile.writeLogAlways(
								"Warning (ODEsolver): Condition of W is "+_W.cond());
					}
					_invW = _W.inverse();
					/*
					 * k1 = invW * ( dYdT(y) + h * d * dFdT )
					 */
					_k1 = _invW.times( _dFdT.times(_h*_d).plus(_dYdT) );
					/*
					 * f1 = dYdT(y + k1*h/2)
					 */
					calc1stDeriv( _k1.times(_h/2).plus(y), _f1 );
					/*
					 * k2 = invW * ( f1 - k1 ) + k1
					 */
					_k2 = _invW.times( _f1.minus(_k1) );
					_k2.plusEquals(_k1);
					/*
					 * ynext = y + h * k2
					 */
					_ynext = _k2.times(_h).plus(y);
					/*
					 * f2 = dYdT(ynext)
					 */
					calc1stDeriv( _ynext, _f2 );
					/*
					 * k3 = invW * ( f2 - e32*(k2-f1) - 2*(k1-y) + h*d*dFdT )
					 * 
					 * First set kaux as the expression inside the brackets,
					 * then multiply by invW on the left.
					 */
					_kaux = _f2.copy();
					_kaux.minusEquals( _k2.minus(_f1).times(_e32) );
					_kaux.minusEquals( _k1.minus(y).times(2) );
					_kaux.plusEquals( _dFdT.times(_h*_d));
					_k3 = _invW.times(_kaux);
					/*
					 * We now use kaux to estimate the error of this step.
					 */
					for (int i = 0; i < _nVar; i++)
						_kaux.set(i, 0, 1/Math.min(y.get(i,0), _ynext.get(i,0)));
					_kaux.arrayTimesEquals(
								_k1.minus(_k2.times(2)).plus(_k3).times(_h/6));
					_error = 0.0;
					for (int i = 0; i < _nVar; i++)
						_error = Math.max(_error, _kaux.get(i,0));
				}
				catch (Exception e)
				{
					LogFile.writeError(e, "Problem in Rosenbrock step");
				}
				/*
				 * The solution is accepted if the weighted error is less than
				 * the relative tolerance rtol. If the step fails, calculate a
				 * new h based on the standard rule for selecting a step size
				 * in numerical integration of initial value problems:
				 * h(n+1) = h(n) * ((rtol / error) ^ power).
				 * 
				 * 90% of this estimated value is then used in the next step to
				 * decrease the probability of further failures.
				 * 
				 * Reference:
				 * GEAR, C. W. 1971. Numerical Initial Value Problems in
				 * Ordinary Differential Equations. Prentice-Hall, Englewood
				 * Cliffs, N.J.
				 */
				_test = Math.pow((rtol/_error), _power);
				if ( _error > rtol )
				{ 
					noFailed = false;
					lastStep = false;
					if ( usingHMin )
						break;
					else if (_EPS * _t > _h * 0.9 * _test)
					{
						usingHMin = true;
						_h = _EPS * _t;
					}
					else
						_h *= 0.9 * _test;
				}
				else
					break;
				LogFile.writeLog("error = "+_error+", rtol = "+rtol+", h = "+_h);
			} // End of `while ( true )`
			/*
			 * If there were no failures compute a new h. We use the same
			 * formula as before to compute a new step, h. But in addition, we
			 * adjust the next time step depending on how stiff the problem is.
			 * If the system is extremely stiff, the increase is limited to
			 * 1.2. Otherwise, the increase is set to a factor of 5.
			 * 
			 * Reference:
			 * Shampine LF. 1982. Implementation of Rosenbrock Methods.
			 * ACM Transactions on Mathematical Software. 8: 93-113.
			 */
			if ( noFailed )
				_h = ( _test < 1.2 ) ? _test : _h * 5;
			/*
			 * The upper limit of hmax still applies.
			 * Update the time.
			 */
			_h = Math.min(_h, hmax);
			_t = _tnext;
			/*
			 * Check no variables have gone negative.
			 * TODO Rob 4June2015: This could be done better.
			 */
			if ( ! _allowNegatives )
				for (int i = 0; i < _nVar; i++)
					if ( _ynext.get(i, 0) < 0.0)
					{
						_ynext.set(i, 0, 0.0);
						LogFile.writeLogAlways(
								"Warning (ODE solver): negative variable! "+i);
					}
			/*
			 * Update the y and the first derivative dYdT.
			 */
			y = _ynext.copy();
			_dYdT = _f2.copy();
		} // End of `while ( ! lastStep )`
		/*
		 * Finally, return the answer.
		 */
		return y;
	}
	
	/**
	 * Update the first derivative of Y, i.e. the rate of change of Y with
	 * respect to time (dYdT = F).
	 * 
	 * @param y
	 * @param target 
	 */
	protected abstract void calc1stDeriv(Matrix y, Matrix target);
	
	/**
	 * Update the second derivative of Y, i.e. the rate of change of F with
	 * respect to time (dFdT).
	 * 
	 * @param y 
	 * @param tdel 
	 * @param target 
	 */
	protected abstract void calc2ndDeriv(Matrix y, Double tdel, Matrix target);
	
	/**
	 * Update the Jacobian matrix, i.e. the rate of change of F with respect to
	 * each of the variables in Y (dFdY).
	 * 
	 * @param y 
	 * @param target 
	 */
	protected abstract void calcJacobian(Matrix y, Matrix target);
	
	/**
	 * 
	 * @param y
	 * @param tdel
	 * @param target
	 */
	protected void numerical2ndDeriv(Matrix y, Double tdel, Matrix target)
	{
		Matrix dYdT = new Matrix(_nVar, 1, 0.0);
		calc1stDeriv(y, dYdT);
		/*
		 * yNext = levels + (tdel * dYdT)
		 */
		Matrix yNext = y.plus(dYdT.times(tdel));
		/*
		 * dFdT = ( dYdT(ynext) - dYdT(y) )/tdel
		 */
		calc1stDeriv(yNext, target);
		target.minusEquals(dYdT);
		target.timesEquals(1/tdel);
	}
}
