/**
 * \package reaction
 * \brief Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS
 * 
 * Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
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

/**
 * \brief Allows creation of a Reaction object whose the reaction rate can be decomposed in several kinetic factors (one factor by solute)
 * 
 * Allows creation of a Reaction object whose the reaction rate can be decomposed in several kinetic factors (one factor by solute)
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class ReactionFactor extends Reaction 
{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Maximum rate at which the reaction may proceed
	 */
	private double            _muMax;
	
	/**
	 * KineticFactor mark-ups define multiplicative factors that affect the reaction rate and decrease the overall rate from muMax 
	 * to something lower. There are several possible types of kinetic factors, and each may take a parameter as well as a 
	 * solute that will set the term's value. This array holds the information in the protocol file for the declared kinetic factor
	 */
	private IsKineticFactor[] _kineticFactor;
	
	/**
	 * A list of the solutes responsible for each kinetic factor. I.e. the solutes which affect this reaction. 
	 * */
	private int[]             _soluteFactor;

	// Temporary variables
	/**
	 * Used to iterate through XML tags that declare this reaction in the protocol file
	 */
	private static int        paramIndex;
	
	/**
	 * Temporary store of values retrieved from the XML protocol file
	 */
	private static double     value;
	
	/**
	 * Marginal rate of reaction matrix
	 */
	private double[]          marginalMu;
	
	/**
	 * Used to calculate diff uptake rate of solute
	 */
	private double[]		  marginalDiffMu;
	
	/**
	 * Holds the unit of the parameter declared in the XML protocol file
	 */
	private StringBuffer      unit;
	
	
	/**
	 * \brief Initialises the reaction by setting relevant parameters to values contained in the simulation protocol file
	 * 
	 * Initialises the reaction by setting relevant parameters to values contained in the simulation protocol file
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param xmlRoot	The XML object containing the definition of one reaction in the protocol file
	 * @see Simulator.createReaction()
	 */
	public void init(Simulator aSim, XMLParser xmlRoot) 
	{

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
			LogFile.writeLog("Error met during ReactionFactor.init(): " + e.getMessage());
		}
	}

	/**
	 * \brief Use the reaction class to fill the parameters fields of the agent
	 * 
	 * Use the reaction class to fill the parameters fields of the agent. Uses information in the agent and protocol file to achieve this
	 * 
	 * @param anAgent	The ActiveAgent which parameters are being populated
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aReactionRoot	The XML object containing the definition of one reaction in the protocol file
	 * @see Simulator.createReaction()
	 */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot) 
	{
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
	 * \brief Return the total number of parameters needed to describe the kinetic of this reaction (muMax included)
	 * 
	 * Return the total number of parameters needed to describe the kinetic of this reaction (muMax included)
	 * 
	 * @return Integer stating the total number of parameters needed to describe the kinetic
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


	/**
	 * \brief Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * @param s	The concentration of solute locally observed
	 * @param mass	Mass of the catalyst (cell...)
	 * @param t	Time
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

	/**
	 * \brief Return the specific reaction rate for a given agent
	 * 
	 * Return the specific reaction rate for a given agent
	 * 
	 * @param anAgent	Agent to use to determine solute concentration and calculate reaction rate
	 * @see ActiveAgent.grow()
	 * @see Episome.computeRate(EpiBac)
	 */
	public void computeSpecificGrowthRate(ActiveAgent anAgent) {

		// Build the array of concentration seen by the agent
		computeSpecificGrowthRate(readConcentrationSeen(anAgent, _soluteList), anAgent);
	}

	/**
	 * \brief Compute specific growth rate in function of concentrations sent Parameters used are those defined for the reaction.
	 * 
	 * Compute specific growth rate in function of concentrations sent Parameters used are those defined for the reaction.
	 * 
	 * @param s	Array of solute concentration
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

	/**
	 * \brief Update the Marginal Mu data matrix
	 * 
	 * Update the Marginal Mu data matrix
	 * 
	 * @param s	Temporary container for solute concentration 
	 */
	public void updateMarginalMu(double[] s) 
	{
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

	/**
	 * \brief Compute the specific growth rate
	 * 
	 * Compute the specific growth rate. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return	The specific growth rate
	 */
	public double computeSpecRate(double[] s){
		double specRate = _muMax;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
			specRate *= marginalMu[iFactor];

		return specRate;
	}

	/**
	 * \brief Compute the marginal difference array
	 * 
	 * Compute the marginal difference array. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return Marginal diff array
	 */
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
	 * \brief Compute specific growth rate in function to concentrations sent
	 * 
	 * Compute specific growth rate in function to concentrations sent
	 * 
	 * @param s	Array of solute concentration
	 * @param anAgent	Parameters used are those defined for THIS agent
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

	/**
	 * \brief Return the marginalDiffMu array
	 * 
	 * Return the marginalDiffMu array. Rob (25/8/2011): added this so that Reaction.applyChemostatReaction() can see marginalDiffMu
	 * 
	 * @return	marginalDiffMu array
	 */
	public double[] getMarginalDiffMu() {
		return marginalDiffMu;
	}

	/**
	 * \brief Calculate the rate of change of each uptake rate with respect to time.
	 * 
	 * Calculate the rate of change of each uptake rate with respect to time. dMUdT = catalyticBiomass*specificGrowthRate*soluteYield.
	 * Returned as a matrix
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param biomass	Total particle mass in the system which catalyses this reaction
	 * @return Matrix containing rate of change of each uptake rate with respect to time
	 */ 
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

	/**
	 * \brief Calculate the rate of change of each uptake rate with respect to each solute. Returned as a matrix
	 * 
	 * Calculate the rate of change of each uptake rate with respect to each solute. Returned as a matrix
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param biomass	Total particle mass in the system which catalyses this reaction
	 * @return Matrix containing rate of change of each uptake rate with respect to each solute
	 */ 
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
	 * \brief Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 *  Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 * @param anAgent	Specific growth rate for this ActiveAgent
	 * @return	The marginal growth rate
	 */
	public double computeMassGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate*anAgent.getParticleMass(_catalystIndex);
	}

	/**
	 * \brief Compute the specific growth rate
	 * 
	 * Compute the specific growth rate
	 * 
	 * @param anAgent	Specific growth rate for this ActiveAgent
	 * @return	The specific growth rate
	 */
	public double computeSpecGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate;
	}
	
}
