/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *______________________________________________________
 * This class allows you to create a Reaction object whose the reaction rate 
 * can be decomposed in several kinetic factor (one factor by solute)
 * 
 */

/**
 * @since January 2007
 * @version 1.0
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */


package simulator.reaction;

import utils.XMLParser;
import utils.UnitConverter;
import utils.LogFile;

import org.jdom.Element;

import Jama.Matrix;

import simulator.Simulator;
import simulator.agent.*;
import simulator.reaction.kinetic.*;

public class ReactionFactor extends Reaction {

	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private double            _muMax;
	private IsKineticFactor[] _kineticFactor;
	/** A list of the solutes responsible for each kinetic factor. I.e. the solutes which affect this reaction. */
	private int[]             _soluteFactor;

	// Temporary variable
	private static int        paramIndex;
	private static double     value;
	private double[]          marginalMu, marginalDiffMu;
	private StringBuffer      unit;

	/* ________________________ CONSTRUCTORS ________________________________ */
	public ReactionFactor() {
	}

	/* ________________ Used during initialisation ______________________ */
	public void init(Simulator aSim, XMLParser xmlRoot) {

		// Call the init of the parent class (populate yield arrays)
		super.init(aSim, xmlRoot);

		// Create the kinetic factors __________________________________________

		// Build the array of different multiplicative limiting expressions
		_kineticFactor = new IsKineticFactor[xmlRoot.getChildren("kineticFactor").size()];
		// one solute factor per kinetic factor
		_soluteFactor = new int[_kineticFactor.length];
		marginalMu = new double[_kineticFactor.length];
		marginalDiffMu = new double[_kineticFactor.length];

		// muMax is the first factor
		unit = new StringBuffer("");
		value = xmlRoot.getParamDbl("muMax", unit);
		_muMax = value*UnitConverter.time(unit.toString());

		int iFactor = 0;
		try {
			// Create and initialise the instance
			for (Element aChild : xmlRoot.getChildren("kineticFactor")) {
				_kineticFactor[iFactor] = (IsKineticFactor) (new XMLParser(aChild))
				.instanceCreator("simulator.reaction.kinetic");
				_kineticFactor[iFactor].init(aChild);
				_soluteFactor[iFactor] = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));
				iFactor++;
			}

			_kineticParam = new double[getTotalParam()];
			_kineticParam[0] = _muMax;

			// Populate the table collecting all kinetic parameters of this
			// reaction term
			paramIndex = 1;
			iFactor = 0;
			for (Element aChild : xmlRoot.getChildren("kineticFactor")) {
				_kineticFactor[iFactor].initFromAgent(aChild, _kineticParam, paramIndex);
				paramIndex += _kineticFactor[iFactor].nParam;
				iFactor++;
			}
		} catch (Exception e) {
			LogFile.writeLog("Error met during ReactionFactor.init()");
		}
	}

	/**
	 * Use the reaction class to fill the parameters fields of the agent
	 */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot) {
		// Call the init of the parent class (populate yield arrays)
		super.initFromAgent(anAgent, aSim, aReactionRoot);

		anAgent.reactionKinetic[reactionIndex] = new double[getTotalParam()];

		// Set muMax
		unit = new StringBuffer("");
		value = aReactionRoot.getParamSuchDbl("kinetic", "muMax", unit);
		double muMax = value*UnitConverter.time(unit.toString());
		anAgent.reactionKinetic[reactionIndex][0] = muMax;

		// Set parameters for each kinetic factor
		paramIndex = 1;
		for (Element aChild : aReactionRoot.getChildren("kineticFactor")) {
			int iSolute = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));
			_kineticFactor[iSolute].initFromAgent(aChild, anAgent.reactionKinetic[reactionIndex],
					paramIndex);
			paramIndex += _kineticFactor[iSolute].nParam;
		}
	}

	/**
	 * @return the total number of parameters needed to describe the kinetic of
	 * this reaction (muMax included)
	 */
	public int getTotalParam() {
		// Sum the number of parameters of each kinetic factor
		int totalParam = 1;
		for (int iFactor = 0; iFactor<_kineticFactor.length; iFactor++) {
			if (_kineticFactor[iFactor]==null) continue;
			totalParam += _kineticFactor[iFactor].nParam;
		}
		return totalParam;
	}

	/* _________________ INTERACTION WITH THE SOLVER_______________________ */

	/**
	 * Update the array of uptake rates and the array of its derivative Based on
	 * parameters sent by the agent
	 * Not called?
	 * @param s
	 * @param mass
	 */
	public void computeUptakeRate(double[] s, ActiveAgent anAgent) {

		// First compute specific rate
		computeSpecificGrowthRate(s, anAgent);

		double mass = anAgent.particleMass[_catalystIndex];

		// Now compute uptake rates
		for (int iSolute : _mySoluteIndex) {
			_uptakeRate[iSolute] = mass*_specRate*anAgent.soluteYield[reactionIndex][iSolute];
		}
		int iSolute;
		for (int i = 0; i<_soluteFactor.length; i++) {
			iSolute = _soluteFactor[i];
			_diffUptakeRate[iSolute] = mass*marginalDiffMu[i]
			                                               *anAgent.soluteYield[reactionIndex][iSolute];
		}
	}

	/**
	 * Update the array of uptake rates and the array of its derivative Based on
	 * default values of parameters Unit is fg.h-1
	 * @param s : the concentration locally observed
	 * @param mass : mass of the catalyst (cell...)
	 */
	public void computeUptakeRate(double[] s, double mass, double t) {

		// First compute specific rate
		computeSpecificGrowthRate(s);

		// Now compute uptake rate
		for (int iSolute : _mySoluteIndex) {
			_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
		}

		int iSolute;
		for (int i = 0; i<_soluteFactor.length; i++) {
			iSolute = _soluteFactor[i];
			if(iSolute!=-1)
				_diffUptakeRate[iSolute] = mass*marginalDiffMu[i]*_soluteYield[iSolute];
		}
	}

	// Rob: created but no longer used
	public void computeUptakeRate(double[] s, double mass, Matrix dFdY) {

		// First compute specific rate
		computeSpecificGrowthRate(s);

		// Now compute uptake rate
		for (int iSolute : _mySoluteIndex) {
			_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
		}

		// Rob (29/9/2011): here we update the chemostat solver's Jacobian matrix directly
		// For each solute factor (a solute which affects this reaction rate), we calculate
		// the effect a change in that solute has on each solute in this reaction's solute
		// index (a list of all solutes which are affected by this reaction).
		int iSolute;
		Matrix soluteYield    = new Matrix (nSolute,1,0);
		Matrix marginaldiffmu = new Matrix (1,nSolute,0);
		// For each affecting solute...
		for (int i = 0; i<_soluteFactor.length; i++){
			iSolute = _soluteFactor[i];
			if(iSolute!=-1){
				marginaldiffmu.set(0, iSolute, marginalDiffMu[i]);
			}
		}
		// For each affected solute...
		for (int i = 0; i<_mySoluteIndex.length; i++){
			iSolute = _mySoluteIndex[i];
			soluteYield.set(iSolute,0,_soluteYield[i]);
		}

		dFdY.plusEquals(soluteYield.times(marginaldiffmu).times(mass));
	}


	/**
	 * Return the specific reaction rate
	 * @see ActiveAgent.grow()
	 * @see Episome.computeRate(EpiBac)
	 */
	public void computeSpecificGrowthRate(ActiveAgent anAgent) {

		// Build the array of concentration seen by the agent
		computeSpecificGrowthRate(readConcentrationSeen(anAgent, _soluteList), anAgent);
	}

	/**
	 * Compute specific growth rate in function of concentrations sent
	 * Parameters used are those defined for the reaction
	 * @param double[] s : array of solute concentration
	 */
	public void computeSpecificGrowthRate(double[] s) {
		_specRate = _muMax;
		int soluteIndex;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1) {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0);
				marginalDiffMu[iFactor] = _muMax*_kineticFactor[iFactor].kineticDiff(0);
			} else {
				marginalMu[iFactor] = _kineticFactor[iFactor]
				                                     .kineticValue(s[_soluteFactor[iFactor]]);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]]);
			}
		}

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor<_soluteFactor.length; jFactor++) {
				if (jFactor!=iFactor) marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
	}

	public void updateMarginalMu(double[] s) {
		int soluteIndex;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1) {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0);
			} else {
				marginalMu[iFactor] = _kineticFactor[iFactor]
				                                     .kineticValue(s[_soluteFactor[iFactor]]);
			}
			//LogFile.writeLog("soluteIndex = "+soluteIndex+", marginalMu = "+marginalMu[iFactor]);
		}
	}

	/** Don't forget to update marginalMu before calling this! */
	public double computeSpecRate(double[] s){
		double specRate = _muMax;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
			specRate *= marginalMu[iFactor];

		return specRate;
	}

	/** Don't forget to update marginalMu before calling this! */
	public double[] computeMarginalDiffMu(double[] s){
		int soluteIndex;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1) {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0);
			} else {
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]]);
			}
			for (int jFactor = 0; jFactor<_soluteFactor.length; jFactor++) {
				if (jFactor!=iFactor) marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}

		return marginalDiffMu;
	}

	/**
	 * Compute specific growth rate in function to concentrations sent
	 * @param double[] s : array of solute concentration
	 * @param anAgent Parameters used are those defined for THIS agent
	 */
	public void computeSpecificGrowthRate(double[] s, ActiveAgent anAgent) {
		double[] kineticParam = anAgent.reactionKinetic[reactionIndex];

		paramIndex = 1;

		// Compute contribution of each limiting solute
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			if (_soluteFactor[iFactor]==-1) { //meaning, if there's no such solute
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0, kineticParam,
						paramIndex);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(0, kineticParam, paramIndex);
			} else {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(
						s[_soluteFactor[iFactor]], kineticParam, paramIndex);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]],
						kineticParam, paramIndex);
			}
			paramIndex += _kineticFactor[iFactor].nParam;
		}

		// First multiplier is muMax
		_specRate = kineticParam[0];

		// Finalise the computation
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor<_soluteFactor.length; jFactor++) {
				if (jFactor!=iFactor) marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
	}

	// Rob (25/8/2011): added this so that Reaction.applyChemostatReaction() can see marginalDiffMu
	public double[] getMarginalDiffMu() {
		return marginalDiffMu;
	}

	/** Rate of change of each uptake rate with respect to time.
	 * dMUdT = catalyticBiomass*specificGrowthRate*soluteYield */
	public Matrix calcdMUdT(Matrix S, double biomass){
		Matrix dMUdT = new Matrix(nSolute, 1, 0);

		try {
			double[] s = new double[nSolute];

			for (int i = 0; i < nSolute; i++) {
				s[i] = S.get(i,0);
				dMUdT.set(i, 0, _soluteYield[i]);
			}

			updateMarginalMu(s);
			_specRate = computeSpecRate(s);

			dMUdT.timesEquals(_specRate*biomass);
		} catch (Exception e) {
			LogFile.writeLog("Error in Reaction.calcdMUdT() : "+e); }
		//LogFile.writeMatrix("dMUdT = ",dMUdT);
		return dMUdT;
	}

	/** Rate of change of each uptake rate with respect to each solute */
	public Matrix calcdMUdS(Matrix S, double biomass){	
		Matrix dMUdY = new Matrix (nSolute, nSolute, 0);

		try {
			double[] s = new double[nSolute];

			for (int i = 0; i < nSolute; i++)
				s[i] = S.get(i,0);

			updateMarginalMu(s);
			marginalDiffMu = computeMarginalDiffMu(s);

			int iSol = 0;
			int jSol = 0;
			// The affecting solute
			for (int iFactor = 0; iFactor < _kineticFactor.length; iFactor++){
				iSol = _soluteFactor[iFactor];
				if (iSol > -1){
					// The affected solute
					for (int jIndex = 0; jIndex < _mySoluteIndex.length; jIndex++){
						jSol = _mySoluteIndex[jIndex];
						dMUdY.set(jSol, iSol, marginalDiffMu[iFactor]*_soluteYield[jSol]);
						//LogFile.writeLog("At "+iSol+", "+jSol+" mDMu = "+marginalDiffMu[iFactor]+" & sY = "+_soluteYield[jSol]);
					}
				}
			}

			dMUdY.timesEquals(biomass);
		} catch (Exception e) {
			LogFile.writeLog("Error in ReactionFactor.calcdMUdS() : "+e); }
		//LogFile.writeMatrix("dMUdY = ", dMUdY);
		return dMUdY;
	}

	/* __________________ Methods called by the agents ___________________ */

	/**
	 * @param anAgent
	 * @return the marginal growth rate (i.e the specific growth rate times the
	 * mass of the particle which is mediating this reaction)
	 */
	public double computeMassGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate*anAgent.getParticleMass(_catalystIndex);
	}

	public double computeSpecGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate;
	}
	
}
