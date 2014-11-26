/**
 * \package reaction
 * \brief Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS
 * 
 * Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.reaction;

import utils.ExtraMath;
import utils.XMLParser;
import utils.UnitConverter;
import utils.LogFile;

import org.jdom.Element;

import Jama.Matrix;
import simulator.Simulator;
import simulator.agent.*;
import simulator.reaction.kinetic.*;

/**
 * \brief Allows creation of a Reaction object whose the reaction rate can be
 * decomposed in several kinetic factor (one factor by solute), but adds a
 * constant factor to the calculation.
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 */
public class ReactionFactorWithConstant extends Reaction 
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Maximum rate at which the reaction may proceed.
	 */
	private Double _muMax;
	
	/**
	 * Constant used to control this reaction.
	 */
	private Double _c;
	
	/**
	 * KineticFactor mark-ups define multiplicative factors that affect the
	 * reaction rate and decrease the overall rate from muMax to something
	 * lower. There are several possible types of kinetic factors, and each may
	 * take a parameter as well as a solute that will set the term's value.
	 * This array holds the information in the protocol file for the declared
	 * kinetic factor.
	 */
	private IsKineticFactor[] _kineticFactor;
	
	/**
	 * A list of the solutes responsible for each kinetic factor. I.e. the
	 * solutes which affect this reaction. 
	 * */
	private int[] _soluteFactor;
	
	// Temporary variables
	private static int iSolute;
	
	/**
	 * Used to iterate through XML tags that declare this reaction in the
	 * protocol file.
	 */
	private static int paramIndex;
	
	/**
	 * Temporary store of values retrieved from the XML protocol file.
	 */
	private static Double value;
	
	/**
	 * Marginal rate of reaction matrix.
	 */
	private Double[] marginalMu;
	
	/**
	 * Used to calculate diff uptake rate of solute.
	 */
	private Double []	 marginalDiffMu;
	
	/**
	 * Holds the unit of the parameter declared in the XML protocol file.
	 */
	private StringBuffer      unit;
	
	/**
	 * \brief Initialises the reaction by setting relevant parameters to values
	 * contained in the simulation protocol file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param xmlRoot The XML object containing the definition of one reaction
	 * in the protocol file.
	 * @see Simulator.createReaction()
	 */
	public void init(Simulator aSim, XMLParser xmlRoot)
	{
		// Call the init of the parent class (populate yield arrays)
		super.init(aSim, xmlRoot);
		// Create the kinetic factors __________________________________________
		// Build the array of different multiplicative limitating expressions
		_kineticFactor = new IsKineticFactor[xmlRoot.getChildrenElements("kineticFactor").size()];
		_soluteFactor = new int[_kineticFactor.length];
		marginalMu = ExtraMath.newDoubleArray(_kineticFactor.length);
		marginalDiffMu = ExtraMath.newDoubleArray(_kineticFactor.length);
		
		// muMax is the first factor
		unit = new StringBuffer("");
		value = xmlRoot.getParamSuchDbl("kinetic", "muMax", unit);
		_muMax = value*UnitConverter.time(unit.toString());
		
		// the constant rate c
		value = xmlRoot.getParamSuchDbl("kinetic", "c", unit);
		_c = value*UnitConverter.time(unit.toString());
		
		int iFactor = 0;
		try
		{
			for (Element aChild : xmlRoot.getChildrenElements("kineticFactor"))
			{
				iSolute = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));
				// Create and initialise the instance
				_kineticFactor[iFactor] = (IsKineticFactor) (new XMLParser(aChild))
									.instanceCreator("simulator.reaction.kinetic.");
				_kineticFactor[iFactor].init(aChild);
				_soluteFactor[iFactor] = iSolute;
				iFactor++;
			}

			_kineticParam = new Double[getTotalParam()];
			_kineticParam[0] = _muMax;
			_kineticParam[1] = _c;

			paramIndex = 2;
			iFactor = 0;
			for (Element aChild : xmlRoot.getChildrenElements("kineticFactor"))
			{
				iSolute = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));
				// Populate the table collecting all kinetic parameters of this
				// reaction term
				_kineticFactor[iFactor].initFromAgent(aChild, _kineticParam, paramIndex);
				paramIndex += _kineticFactor[iFactor].nParam;
				iFactor++;
			}
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "ReactionFactor.init()");
		}
	}
	
	/**
	 * \brief Use the reaction class to fill the parameters fields of the agent.
	 * 
	 * Uses information in the agent and protocol file to achieve this.
	 * 
	 * @param anAgent	The ActiveAgent which parameters are being populated.
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aReactionRoot	The XML object containing the definition of one
	 * reaction in the protocol file.
	 * @see Simulator.createReaction()
	 */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot)
	{
		// Call the init of the parent class (populate yield arrays)
		super.initFromAgent(anAgent, aSim, aReactionRoot);
		anAgent.reactionKinetic[reactionIndex] = ExtraMath.newDoubleArray(getTotalParam());
		
		// Set muMax
		unit = new StringBuffer("");
		value = aReactionRoot.getParamSuchDbl("kinetic", "muMax", unit);
		Double muMax = value*UnitConverter.time(unit.toString());
		anAgent.reactionKinetic[reactionIndex][0] = muMax;
		
		// Set c
		value = aReactionRoot.getParamSuchDbl("kinetic", "c", unit);
		Double c = value*UnitConverter.time(unit.toString());
		anAgent.reactionKinetic[reactionIndex][1] = c;
		
		// Set parameters for each kinetic factor (muMax and c are the first 2)
		paramIndex = 2;
		for (Element aChild : aReactionRoot.getChildrenElements("kineticFactor"))
		{
			iSolute = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));
			_kineticFactor[iSolute].initFromAgent(aChild,
						anAgent.reactionKinetic[reactionIndex], paramIndex);
			paramIndex += _kineticFactor[iSolute].nParam;
		}
	}
	
	/**
	 * \brief Return the total number of parameters needed to describe the
	 * kinetic of this reaction (muMax included).
	 * 
	 * @return Integer stating the total number of parameters needed to
	 * describe the kinetic.
	 */
	public int getTotalParam()
	{
		// Sum the number of parameters of each kinetic factor
		int totalParam = 2;
		for (int iFactor = 0; iFactor < _kineticFactor.length; iFactor++)
		{
			if (_kineticFactor[iFactor] == null)
				continue;
			totalParam += _kineticFactor[iFactor].nParam;
		}
		return totalParam;
	}
	
	/**
	 * \brief Update the array of uptake rates and the array of its derivative.
	 * 
	 * Based on default values of parameters. Unit is fg.h-1
	 * 
	 * @param s	The concentration of solute locally observed
	 * @param mass	Mass of the catalyst (cell...)
	 * @param tdel	Time
	 */
	public void computeUptakeRate(Double[] s, Double mass, Double tdel)
	{
		// First compute specific rate
		computeSpecificGrowthRate(s);
		
		if(Simulator.isChemostat)
		{
			for (int iSolute : _mySoluteIndex)
				_uptakeRate[iSolute] = (tdel*mass*Dil) + (mass *_specRate*_soluteYield[iSolute]);
			int iSolute;
			for (int i = 0; i<_soluteFactor.length; i++)
			{
				iSolute = _soluteFactor[i];
				if(iSolute!=-1)
					_diffUptakeRate[iSolute] = (tdel*mass*Dil) + (mass*marginalDiffMu[i]*_soluteYield[iSolute]);
			}

		}
		else
		{
			for (int i = 0; i < _mySoluteIndex.length; i++)
			{
				iSolute = _mySoluteIndex[i];
				_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
			}

			for (int i = 0; i < _soluteFactor.length; i++)
			{
				iSolute = _soluteFactor[i];
				_diffUptakeRate[iSolute] = mass*marginalDiffMu[i]*_soluteYield[iSolute];
			}
		}
	}
	
	/**
	 * \brief Return the specific reaction rate for a given agent.
	 * 
	 * @param anAgent	Agent to use to determine solute concentration and
	 * calculate reaction rate.
	 * @see ActiveAgent.grow()
	 * @see Episome.computeRate(EpiBac)
	 */
	public void computeSpecificGrowthRate(ActiveAgent anAgent)
	{
		computeSpecificGrowthRate(readConcentrationSeen(anAgent, _soluteList),anAgent);
	}
	
	/**
	 * \brief Compute specific growth rate in function of concentrations sent.
	 * 
	 * Parameters used are those defined for the reaction.
	 * 
	 * @param s	Array of solute concentration
	 */
	public void computeSpecificGrowthRate(Double[] s)
	{
		_specRate = _muMax;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
		{
			iSolute = _soluteFactor[iFactor];
			marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(s[iSolute]);
			marginalDiffMu[iFactor] = _muMax*_kineticFactor[iFactor].kineticDiff(s[iSolute]);
		}
		for (int iFactor = 0; iFactor < _soluteFactor.length; iFactor++)
		{
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor < _soluteFactor.length; jFactor++)
			{
				if (jFactor != iFactor)
					marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
		_specRate += _c; 
	}
	
	/**
	 * \brief Compute specific growth rate in function to concentrations sent.
	 * 
	 * @param s	Array of solute concentration.
	 * @param anAgent	Parameters used are those defined for THIS agent.
	 */
	public void computeSpecificGrowthRate(Double[] s, ActiveAgent anAgent)
	{
		Double[] kineticParam = anAgent.reactionKinetic[reactionIndex];
		// First multiplier is muMax
		_specRate = kineticParam[0];
		paramIndex = 2;
		// Compute contribution of each limiting solute
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
		{
			iSolute = _soluteFactor[iFactor];
			marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(
										s[iSolute], kineticParam, paramIndex);
			marginalDiffMu[iFactor] = _muMax *
							_kineticFactor[iFactor].kineticDiff(s[iSolute],
													kineticParam, paramIndex);
			paramIndex += _kineticFactor[iFactor].nParam;
		}
		
		// Finalize the computation
		for (int iFactor = 0; iFactor < _soluteFactor.length; iFactor++)
		{
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor < _soluteFactor.length; jFactor++)
			{
				if (jFactor!=iFactor)
					marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
		// add constant rate
		_specRate += kineticParam[1]; 
	}
	
	/**
	 * \brief Calculate the rate of change of each uptake rate with respect to
	 * each solute.
	 * 
	 * @param S	Temporary container for solute concentration.
	 * @param biomass	Total particle mass in the system which catalyses this
	 * reaction.
	 * @return Matrix containing rate of change of each uptake rate with respect
	 * to each solute.
	 */ 
	public Matrix calcdMUdS(Matrix S, Double biomass)
	{
		Matrix dMUdY = new Matrix (nSolute, nSolute, 0.0);
		try
		{
			// No need to initialise this with ExtraMath.newDoubleArray(nI)
			Double[] s = new Double[nSolute];
			for (int i = 0; i < nSolute; i++)
				s[i] = S.get(i, 0);
			
			updateMarginalMu(s);
			marginalDiffMu = computeMarginalDiffMu(s);
			
			int iSol, jSol;
			// The affecting solute
			for (int iFactor = 0; iFactor < _kineticFactor.length; iFactor++)
			{
				iSol = _soluteFactor[iFactor];
				if (iSol > -1)
				{
					// The affected solute
					for (int jIndex = 0; jIndex < _mySoluteIndex.length; jIndex++)
					{
						jSol = _mySoluteIndex[jIndex];
						dMUdY.set(jSol, iSol, marginalDiffMu[iFactor]*_soluteYield[jSol]);
					}
				}
			}
			dMUdY.timesEquals(biomass);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "ReactionFactorWithConstant.calcdMUdS()");
		}
		return dMUdY;
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
	public Matrix calcdMUdT(Matrix S, Double biomass)
	{
		Matrix dMUdT = new Matrix(nSolute, 1, 0.0);
		try
		{
			// No need to initialise this with ExtraMath.newDoubleArray(nI)
			Double[] s = new Double[nSolute];
			for (int i = 0; i < nSolute; i++)
			{
				s[i] = S.get(i, 0);
				dMUdT.set(i, 0, _soluteYield[i]);
			}
			
			updateMarginalMu(s);
			_specRate = computeSpecRate(s);
			
			dMUdT.timesEquals(_specRate*biomass);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "Reaction.calcdMUdT()");
		}
		return dMUdT;
	}

	/**
	 * \brief Compute the marginal difference array
	 * 
	 * Compute the marginal difference array. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return Marginal diff array
	 */
	public Double[] computeMarginalDiffMu(Double[] s) {
		int soluteIndex;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1) {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0);
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
	 * \brief Compute the specific growth rate
	 * 
	 * Compute the specific growth rate. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return	The specific growth rate
	 */
	public Double computeSpecRate(Double[] s) {
		Double specRate = _muMax;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
			specRate *= marginalMu[iFactor];

		return specRate;
	}

	/* __________________ Methods called by the agents ___________________ */

	/**
	 * \brief Update the Marginal Mu data matrix
	 * 
	 * Update the Marginal Mu data matrix
	 * 
	 * @param s	Temporary container for solute concentration 
	 */
	public void updateMarginalMu(Double[] s) {
		int soluteIndex;

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1) {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0);
			} else {
				marginalMu[iFactor] = _kineticFactor[iFactor]
				                                     .kineticValue(s[_soluteFactor[iFactor]]);
			}
			//LogFile.writeLog("soluteIndex = "+soluteIndex+", marginalMu = "+marginalMu[iFactor]);
		}

	}

	/**
	 * \brief Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 *  Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 * @param anAgent	Specific growth rate for this ActiveAgent
	 * @return	The marginal growth rate
	 */
	public Double computeMassGrowthRate(ActiveAgent anAgent) {
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
	public Double computeSpecGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate;
	}
}