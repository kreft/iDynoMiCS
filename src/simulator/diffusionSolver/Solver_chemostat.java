/**
 * \package diffusionSolver
 * \brief Package of classes used to capture the diffusion solvers that can be defined in the protocol file
 * 
 * Package of classes used to capture the diffusion solvers that can be defined in the protocol file. Solvers are used to compute 
 * the steady-state solute profile within the computational domains. This package is  part of iDynoMiCS v1.2, governed by the CeCILL 
 * license under French law and abides by the rules of distribution of free software. You can use, modify and/ or redistribute iDynoMiCS 
 * under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.diffusionSolver;

import java.util.ArrayList;
import Jama.Matrix;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.diffusionSolver.multigrid.MultigridSolute;
import simulator.geometry.Bulk;
import simulator.geometry.Domain;
import simulator.geometry.boundaryConditions.AllBC;
import utils.LogFile;
import utils.XMLParser;
import idyno.Idynomics;
import idyno.SimTimer;

public class Solver_chemostat extends DiffusionSolver 
{
	/** 
	 * Number of solutes in simulation 
	 * */
	protected int nSolute;
	
	/** 
	 * Number of reactions in simulation 
	 * */
	protected int nReaction;

	/**
	 * The solute grid of solute concentrations. As we are in a chemostat this has only one grid position (0,0,0). Note that y will be 
	 * the vector containing the actual concentrations during odeSolver()... here allSolute will be a temporary variable used for 
	 * calculating potential reaction rates, etc, in calcReacAndDiffReacRates().
	 */
	protected  SoluteGrid[]      	allSolute;
	
	/**
	 * The solute grid of reaction rates for each solute. 
	 * */
	protected  SoluteGrid[] 		allReac;
	
	/**
	 * Computational domain that this solver is associated with
	 */
	protected Domain                _domain;
	
	/**
	 * Array of reactive biomasses
	 */
	protected MultigridSolute []    _reactiveBiomass;
	
	/**
	 * Boolean stating whether the bulk concentration remains constant
	 */
	protected Boolean[] isConstSol;

	/**
	 * Constant used by the ODE solvers
	 */
	protected double d    = 1.0 / (2.0 + Math.sqrt(2.0));
	
	/**
	 * Constant used by the ODE solvers
	 */
	protected double e32  = 6.0 + Math.sqrt(2.0);
	
	/**
	 * Constant used by the ODE solvers
	 */
	protected double power = 1.0/3.0;
	
	/**
	 * Numerical accuracy for EPS (error per step) 
	 */
	protected double sqrtE = Math.sqrt(2.22e-16);
	
	/**
	 * The smallest positive floating-point number such that 1.0+EPS > 1.0 
	 */
	protected double EPS = 2.22e-16;
	
	/**
	 * The chemostat dilution rate. This will be set when init() calls setDilutionAndY0() 
	 */
	protected double Dilution; 

	/**
	 * Maximum internal step of the solver. Used in the ODE solver
	 */
	protected double hmax;
	
	/**
	 * Relative tolerance of the calculated error. Used in the ODE solver
	 */
	protected double rtol;

	/** 
	 * 1D array of arraylists containing the list of reactions in which each of the solutes participates. 
	 */
	protected ArrayList<Integer> solReacInd = new ArrayList<Integer>();
	
	/** 
	 * For each reaction (first index) this stores the yield of each solute (second index). Used in calculating the Jacobian. The 
	 * Jacobian matrix: on each row an "F" (the rate of change of a substrate concentration) is differentiated with respect to each 
	 * of the substrate concentrations (arranged on the columns).
	 */ 
	protected double[][] soluteYield;
	
	


	/**
	 * \brief Initialise the chemostat: the solute and reaction list, and all variables used by the ODE solvers
	 * 
	 * Initialise the chemostat: the solute and reaction list, and all variables used by the ODE solvers
	 * 
	 * @param aSimulator	The current simulation object
	 * @param xmlRoot	Root element of XML tags containing chemostat related parameters
	 */
	public void init(Simulator aSimulator, XMLParser xmlRoot) 
	{

		super.init(aSimulator, xmlRoot);

		// _soluteList and _reactions are in DiffusionSolver
		nSolute = _soluteList.length;
		nReaction = _reactions.size();
		//LogFile.writeLog(" nSolute = "+nSolute+", nReaction = "+nReaction);

		allSolute = new SoluteGrid[nSolute];
		allReac = new SoluteGrid[nSolute];

		// We read in the domain from the protocol file
		_domain = aSimulator.world.getDomain(xmlRoot.getAttribute("domain"));

		// Initialise variables used by the ODE solvers
		hmax = xmlRoot.getParamDbl("hmax");
		rtol = xmlRoot.getParamDbl("rtol");
		setDilutionAndY0();

		// The soluteYield (constant) and allDiffReac (variable) will be used for 
		// calculating the Jacobian matrix
		soluteYield = new double [nReaction][nSolute];
		try{
			for (int i = 0; i < nReaction; i++)
				soluteYield[i] = _reactions.get(i).getSoluteYield();
		}catch(Exception e){LogFile.writeLogAlways("Error creating solute fields in Solver_chemostat.init()");}
		// Initialize array of reactive biomasses
		_reactiveBiomass = new MultigridSolute[nReaction];
		try{
			for (int i = 0; i<nReaction; i++) {
				_reactiveBiomass[i] = new MultigridSolute(_soluteList[0], _reactions.get(i).reactionName);
				_reactiveBiomass[i]._conc[0].setAllValueAt(0);
				if (!Idynomics.quietMode)
					System.out.println("biomass conc is ----->>>   " + _reactiveBiomass[i]._conc[0].grid[0][0][0]);
			}
		}catch(Exception e){LogFile.writeLogAlways("Error creating reaction fields in Solver_chemostat.init()");}
	}

	/**
	 * \brief Set the dilution rate and the initial substrate concentration to those specified in the bulk compartment
	 * 
	 * Set the dilution rate and the initial substrate concentration to those specified in the bulk compartment 
	 * [Rob: one assumption here is that the dilution rate is constant].
	 */
	public void setDilutionAndY0() {
		try {
			for (AllBC aBC : _domain.getAllBoundaries()){
				if (aBC.hasBulk()){
					Bulk aBulk = aBC.getBulk();
					if(aBulk.getName().equals("chemostat")){
						Dilution = aBulk._D;
						for (int i = 0; i < nSolute; i++) {
							allSolute[i] = mySim.soluteList[i];
							allSolute[i].setAllValueAt(aBulk._bulkValue[i]);
							//LogFile.writeLog(" allSolute["+i+"] = "+allSolute[i].grid[0][0][0]);
							allReac[i] = new SoluteGrid (1,1,1, _domain._resolution,
									mySim.soluteList[i].gridName, mySim.soluteList[i].getDomain());
							allReac[i].setAllValueAt(0);
						}
					}
					isConstSol = aBulk._isConstant;
				}
				//LogFile.writeLog("S = "+allSolute[0].grid[0][0][0]+" X = "+allSolute[1].grid[0][0][0]);
			}
		} catch (Exception e) {
			LogFile.writeLogAlways("Error in Solver_chemostat.updateReacRateAndDiffRate() : " + e);}
	}

	/**
	 * \brief Initialise the concentration fields within the chemostat, getting concentration of any catalysts and storing this on a grid
	 * 
	 * Initialise the concentration fields within the chemostat, getting concentration of any catalysts and storing this on a grid
	 */
	public void initializeConcentrationFields() {
		try {
			//reset biomass concentration in the grid
			for (int i = 0; i<nReaction; i++) {
				_reactiveBiomass[i]._conc[0].setAllValueAt(0);
			}

			// Get the catalyst (biomass and other particulates) CONCENTRATION
			for (int i = 0; i<nReaction; i++) {
				_reactions.get(i).fitAgentMassOnGrid(_reactiveBiomass[i].getGrid());
				if (!Idynomics.quietMode)
					System.out.println("biomass conc is ----->>>   " + _reactiveBiomass[i]._conc[0].grid[0][0][0]);
			}
		} catch (Exception e) {
			LogFile.writeLogAlways("Error in Solver_chemostat.initializeConcentrationFields() : " + e);}
	}

	/**
	 * \brief Use the ODE solver to solve the diffusion reaction and update the bulk 
	 * 
	 * Use the ODE solver to solve the diffusion reaction and update the bulk
	 * 
	 */
	public void solveDiffusionReaction() 
	{
		//LogFile.writeLog("S = "+allSolute[0].grid[0][0][0]+" X = "+allSolute[1].grid[0][0][0]);
		odeSolver(SimTimer.getCurrentTime(), rtol, hmax);
		updateBulk();
	}

	/**
	 * \brief ODE solver for calculating the diffusion reactions
	 * 
	 * ODE solver for calculating the diffusion reactions
	 * 
	 * @param t0	Simulation time
	 * @param rtol	Relative tolerance of the calculated error
	 * @param hmax	Maximum internal step of the solver
	 */
	public void odeSolver(double t0, double rtol, double hmax) {

		Matrix y       = new Matrix(nSolute,1,0);
		Matrix ynext   = new Matrix(nSolute,1,0);
		double t, tnext, tfinal;
		double h, error;
		Matrix f1      = new Matrix(nSolute,1,0);
		Matrix f2      = new Matrix(nSolute,1,0);
		Matrix W       = new Matrix (nSolute, nSolute, 0);
		Matrix invW    = new Matrix (nSolute, nSolute, 0);
		Matrix k1      = new Matrix (nSolute, 1, 0);
		Matrix k2      = new Matrix (nSolute, 1, 0);
		Matrix k3      = new Matrix (nSolute, 1, 0);
		Matrix kaux    = new Matrix (nSolute, 1, 0);
		Matrix dFdY    = new Matrix (nSolute, nSolute, 0);
		Matrix dYdT    = new Matrix (nSolute, 1, 0);
		Matrix dFdT    = new Matrix (nSolute, 1, 0);
		Matrix sInflow = new Matrix (nSolute, 1, 0); // This is counted as a variable as solute may be pulsed


		// We reset t and tnext to zero and find how long the current time-step is
		t = 0;
		tnext = 0;
		h = 0;
		tfinal =  SimTimer.getCurrentTimeStep();
		// The error estimate is set back to zero for each global time-step
		error = 0;

		for (int iSol = 0; iSol<nSolute; iSol++) {
			y.set(iSol, 0,  allSolute[iSol].grid[0][0][0]);
		}
		
		// Check if the Sinflow has changed (solutes may be pulsed)
		sInflow = updateSInflow(sInflow);

		dYdT = calcdYdT(y, sInflow, dYdT);


		// Control statement in case hmax > global time-step
		while (hmax > tfinal) {
			hmax *= 0.5;
			rtol *= 0.5;
		}

		// hmax is our first attempted step - if it is too large then it will be reduced
		h = hmax;



		//>>>>>>>>>>---> LOOP for h-steps <---<<<<<<<<<<//
		boolean lastStep = false;

		while (!lastStep) {

			// If a step is successful then h will be increased:
			// this stops it from growing too large
			h = Math.min(hmax,h);

			// If the time step is within 5% of tfinal, take h = tfinal-t;
			if (t + (1.05 * h) >= tfinal){
				h = tfinal - t;
				lastStep = true;
			}

			// tdel is a mini time-step used for calculating the value of dFdT 
			double tdel = sqrtE*(t+h);
			dFdT = calcdFdT(y, dFdT, sInflow, tdel);
			//LogFile.writeMatrix("t = "+t+", dFdT = ", dFdT);

			// Update the uptake rates and diffuptake rates and then the dFdY Jacobian matrix
			//updateDerivatives(y, dFdY);
			dFdY = calcJacobian(y,dFdY);
			//LogFile.writeMatrix("dFdY = ",dFdY);

			//>>>---> LOOP for trying out this h <---<<<//
			// 
			boolean noFailed = true;
			boolean usingHMin = false;
			while (true) {

				//LogFile.writeLog("h = "+h);
				if (lastStep){
					tnext = tfinal;
				} else {
					tnext = t + h;
				}

				// The Rosenbrock method
				//LogFile.writeLog("Starting Rosenbrock. dFdY(0,0) = "+dFdY.get(0, 0));
				try {
					// W = I - h * d * dFdY
					W = dFdY.times(-h*d);
					//LogFile.writeLog("W(0,0) = "+W.get(0, 0));
					W.plusEquals(Matrix.identity(nSolute,nSolute));
					if (W.cond() > 10)
						LogFile.writeLogAlways("Warning (ODEsolver): Condition of W is "+W.cond());
					invW = W.inverse();

					// k1 = invW * ( dYdT(y) + h * d * dFdT )
					k1.timesEquals(0);
					k1 = invW.times( dFdT.times(h*d).plus(dYdT) );
					//LogFile.writeLog("k1 = "+k1.get(0,0));

					// f1 = dYdT(y + k1*h/2)
					f1.timesEquals(0);
					f1 = calcdYdT( k1.times(h/2).plus(y) , sInflow, f1);
					//LogFile.writeLog("f1 = "+f1.get(0,0));

					// k2 = invW * ( f1 - k1 ) + k1
					k2.timesEquals(0);
					k2 = invW.times( f1.minus(k1) );
					k2.plusEquals(k1);
					//LogFile.writeLog("k2 = "+k2.get(0,0));

					// ynext = y + h * k2
					ynext.timesEquals(0);
					ynext = k2.times(h).plus(y);
					//LogFile.writeLog("ynext = "+ynext.get(0,0));

					// f2 = dYdT(ynext)
					f2.timesEquals(0);
					f2 = calcdYdT(ynext, sInflow, f2);
					//LogFile.writeLog("f2 = "+f2.get(0,0));

					// k3 = invW * ( f2 - e32*(k2-f1) - 2*(k1-y) + h*d*dFdT )
					// Setting kaux as expression inside brackets
					k3.timesEquals(0);
					kaux.timesEquals(0);
					kaux = f2.copy();
					kaux.minusEquals( k2.minus(f1).times(e32) );
					kaux.minusEquals( k1.minus(y).times(2) );
					kaux.plusEquals( dFdT.times(h*d));
					k3 = invW.times(kaux);
					//LogFile.writeLog("k3 = "+k3.get(0,0));
					// error = (h/6) * (k1 - 2*k2 + k3)/y
					kaux.timesEquals(0);
					// We now use kaux to estimate the error of this step

					for (int i = 0; i < nSolute; i++)
						kaux.set(i,0, 1/Math.min(y.get(i,0),ynext.get(i,0)));
					//LogFile.writeLog("y "+y.get(0, 0)+", ynext "+ynext.get(0,0)+", kaux "+kaux.get(0,0));
					kaux.arrayTimesEquals( k1.minus( k2.times(2) ).plus(k3).times(h/6) );
					//kaux = k1.minus( k2.times(2) ).plus(k3).times(h/6);
					// We now calculate the error
					error = 0;
					for (int i = 0; i < nSolute; i++)
						error = Math.max(error, kaux.get(i,0));
					//LogFile.writeLog("error = "+error);

				} catch (Exception e) {
					LogFile.writeLogAlways("Error in Solver_chemostat.rosenbrock() : " + e);}


				// The solution is accepted if the weighted error is less than the relative tolerance rtol.
				if(error > rtol) { 

					noFailed = false;
					lastStep = false;		

					// If the step fails, calculate a new h based on the standard rule for 
					// selecting a step size in numerical integration of initial value problems: 
					// hn+1 = (rtol / error) ^ (1/order of method,in our case is 3) * hn;
					// 90% of this estimated value is then used in the next step to
					// decrease the probability of further failures.

					// Reference: GEAR, C. W. 1971. Numerical Initial Value Problems in 
					// Ordinary Differential Equations. Prentice-Hall, Englewood Cliffs, N.J.

					if (usingHMin){
						break;
					} else {
						if (EPS*t > h * 0.9*(Math.pow((rtol/error), power)))
							usingHMin = true;
						h =  Math.max(EPS*t, h * 0.9*(Math.pow((rtol/error), power)));
					}
					// Note that EPS*t = hmin (variable previously used)


				}else{
					//LogFile.writeLog("tnext = "+t+", ynext = "+ynext.get(1,0)+", dYdT = "+dYdT.get(1,0));
					break;
				}// End of if(error > rtol)
				LogFile.writeLog("error = "+error+", rtol = "+rtol+", h = "+h);
			}// End of while(true)

			//If there were no failures compute a new h.
			if(noFailed) {
				// We use the same formula as before to compute a new step, h. But in addition,
				// we adjust the next time step depending on how stiff the problem is.

				// Reference: Shampine LF. 1982. Implementation of Rosenbrock Methods. 
				// ACM Transactions on Mathematical Software. 8: 93-113.
				double test = Math.pow((rtol/error), power);

				if(test<1.2) {
					// If the system is extremely stiff, the increase is limited to 1.2
					h *= test;
				} else {
					// Otherwise, the increase is set to a factor of 5
					h *= 5;
				}


			}// End of if(noFailed)

			// Update the time and the substrate concentration
			t = tnext;

			for (int iSol = 0; iSol < nSolute; iSol++){
				// For each solute we check first whether its concentration should remain constant
				if (!isConstSol[iSol]){
					// We then see if its concentration has gone negative (this is a problem!)
					if (y.get(iSol,0) < 0){
						y.set(iSol, 0, 0);
						LogFile.writeLogAlways("Warning! Solute has gone negative");
					} else {
						y.set(iSol, 0, ynext.get(iSol, 0));
					}
				}
			}// End of for(int iSol=0; iSol<nSolute; iSol++)

			dYdT = f2.copy();
			//LogFile.writeMatrix("y = ", y);
			//LogFile.writeMatrix("dYdT = ", dYdT);

		}// End of while(!lastStep)

		//LogFile.writeMatrix("y = ", y);
		//LogFile.writeMatrix("dYdT = ", dYdT);
		//LogFile.writeMatrix("dFdT = ", dFdT);
		//LogFile.writeMatrix("dFdY = ", dFdY);
		for (int iSol = 0; iSol < nSolute; iSol++)
			// Assuming all is well, we then update allSolute to the appropriate value
			allSolute[iSol].setAllValueAt(y.get(iSol,0));
	}

	/**
	 * \brief Calculates the derivative of substrate concentration (a function parameter) with respect to time. 
	 * 
	 * Calculates the derivative of substrate concentration (a function parameter) with respect to time. As part of this process the 
	 * allReac and allDiffReac values are set according to the (possibly hypothetical) values at substrate concentration S - you'll 
	 * probably want to rerun calcReacRateAndDiffRate(y); afterwards to get allReac and allDiffReac to the correct values.
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param sInflow	Matrix of solute flow information (if solutes are pulsed)
	 * @param dYdT	Matrix to store the derivatives of substrate concentration
	 * @return Matrix of the derivatives of substrate concentration
	 */
	public Matrix calcdYdT(Matrix S, Matrix sInflow, Matrix dYdT) 
	{
		try 
		{
			dYdT = sInflow.minus(S).times(Dilution);
			//LogFile.writeMatrix("D(Sin - S) = ", dYdT);

			for (int iReac = 0; iReac < nReaction; iReac++){
				//LogFile.writeLog("Reaction "+iReac);
				dYdT.plusEquals(_reactions.get(iReac).calcdMUdT(S, _reactiveBiomass[iReac]._conc[0].grid[0][0][0]));
			}
		} catch (Exception e) {
			LogFile.writeLogAlways("Error in Solver_chemostat.calcdYdT() : " + e);}
		//LogFile.writeMatrix("calcdYdT makes dYdT to be", dYdT);
		return dYdT;
	}

	/** \brief Provides an estimate of how F (the derivative of S w.r.t. time) is changing (i.e. the second derivative).
	 * 
	 * Provides an estimate of how F (the derivative of S w.r.t. time) is changing (i.e. the second derivative). Note that calcdYdT 
	 * is called with the hypothetical S-value of y + dYdT * tdel, so we'll need to recall calcReacRateAndDiffRate(S) with S as just y.
	 * 	
	 * @param S	Temporary container for solute concentration
	 * @param dFdT	Estimation of how the derivative of solutes with regard to time is changing
	 * @param sInflow	Matrix of solute flow information (if solutes are pulsed)
	 * @param tdel	Mini time-step used for calculating the value of dFdT 
	 */
	public Matrix calcdFdT(Matrix S, Matrix dFdT, Matrix sInflow, double tdel) 
	{
		Matrix Snext = new Matrix (nSolute, 1, 0);
		Matrix dYdT = new Matrix (nSolute, 1, 0);

		// ynext = y + tdel * dYdT
		dYdT = calcdYdT(S, sInflow, dYdT);
		Snext = S.plus(dYdT.times(tdel));

		// dFdT = ( dYdT(ynext) - dYdT(y) )/tdel
		dFdT = calcdYdT(Snext, sInflow, dFdT);
		dFdT.minusEquals(dYdT);
		dFdT.timesEquals(1/tdel);

		// We recalculate dYdT to make sure that marginalMu values are correct
		dYdT = calcdYdT(S, sInflow, dYdT);

		return dFdT;
	}

	/**
	 * \brief Calculate the Jacobian matrix - the partial derivatives of the rate of change of each substrate
	 * 
	 * The Jacobian matrix dFdY gives the partial derivatives of the rate of change of each substrate (F = dYdT, or equivalently, dSdT) 
	 * with respect to each substrate. Rows correspond to substrates being differentiated (first by time to give F, then by other 
	 * substrate concentrations to give the elements of this matrix) and columns to the substrate with respect to which we are 
	 * differentiating.
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param dFdY	Jacobian matrix being calculated
	 * @return Matrix containing partial derivatives of rate of change of each substrate
	 */
	public Matrix calcJacobian(Matrix S, Matrix dFdY) {

		double biomass = 0;
		dFdY = Matrix.identity(nSolute,nSolute).times(-Dilution);

		for (int iReac = 0; iReac < nReaction; iReac++){
			//LogFile.writeLog("Reaction "+iReac);
			// biomass is the total particle mass in the system which catalyses this reaction
			biomass = _reactiveBiomass[iReac]._conc[0].grid[0][0][0];
			dFdY.plusEquals(_reactions.get(iReac).calcdMUdS(S, biomass));
		}

		return dFdY;
	}

	/**
	 * \brief Find the connected bulks and update their concentrations
	 * 
	 * Find the connected bulks and update their concentrations
	 */
	public void updateBulk() 
	{
		try {
			for (AllBC aBC : myDomain.getAllBoundaries())
				if (aBC.hasBulk()) aBC.updateBulk(allSolute, allReac, 0);
		} catch (Exception e) {
			LogFile.writeLog("Error in Solver_chemostat.updateBulk() : " + e);}
	}

	/**
	 * \brief Check if the Sinflow has changed (solutes may be pulsed)
	 * 
	 * Check if the Sinflow has changed (solutes may be pulsed)
	 * 
	 * @param sInflow	Matrix of solutes in flow
	 * @return	Updated matrix showing solute flow levels
	 */
	public Matrix updateSInflow(Matrix sInflow) {
		try {

			for (AllBC aBC : _domain.getAllBoundaries()){
				if (aBC.hasBulk()){
					Bulk aBulk = aBC.getBulk();
					if(aBulk.getName().equals("chemostat")){
						for (int i = 0; i < nSolute; i++)
							sInflow.set(i, 0, aBulk._sIn[i]);
					}
				}	
			}

		} catch (Exception e) {
			LogFile.writeLog("Error in Solver_chemostat.updateSInflow() : " + e);}

		return sInflow;

	}

}
