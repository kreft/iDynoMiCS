/**
 * \package reaction
 * \brief Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS
 * 
 * Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.reaction;

import idyno.SimTimer;

import java.util.ArrayList;

import org.jdom.Element;

import Jama.Matrix;
import simulator.Simulator;
import simulator.agent.ActiveAgent;
import simulator.agent.zoo.MultiEpiBac;
import simulator.agent.zoo.MultiEpisomeParam;
import simulator.reaction.kinetic.IsKineticFactor;
import utils.ExtraMath;
import utils.LogFile;
import utils.UnitConverter;
import utils.XMLParser;

/**
 * \brief Captures the fitness cost associated with the carriage of a plasmid
 * 
 * This class deals with the fitness cost associated with the carriage of a plasmid. A cell carrying a plasmid will have its biomass 
 * yield reduced due to the metabolic burden conferred by the plasmid, which indirectly impacts its growth rate output. Here we assume 
 * that plasmid cost decreases exponentially with the time spent in the host, until it reaches a minimum value. The assumption behind 
 * this approach is that plasmid and host coevolve. We also assume that the individuals from a species don't mutate, thus a plasmid 
 * that has evolved in a host of a certain species if transferred to another individual of the same species, will have the "evolved" 
 * cost for its new recipient.
 * 
 * @author SoniaMartins
 *
 */
public class ReactionGrowthFitness extends Reaction{

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
	private Double[] marginalDiffMu;
	
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
	@Override
	public void init(Simulator aSim, XMLParser xmlRoot) {

		// Call the init of the parent class (populate yield arrays)
		super.init(aSim, xmlRoot);

		// Create the kinetic factors __________________________________________

		// Build the array of different multiplicative limiting expressions
		_kineticFactor = new IsKineticFactor[xmlRoot.getChildrenElements("kineticFactor").size()];
		_soluteFactor = new int[_kineticFactor.length];
		marginalMu = new double[_kineticFactor.length];
		marginalDiffMu = ExtraMath.newDoubleArray(_kineticFactor.length);

		// muMax is the first factor
		unit = new StringBuffer("");
		value = xmlRoot.getParamDbl("muMax", unit);
		_muMax = value*UnitConverter.time(unit.toString());

		int iFactor = 0;
		int iSolute;
		try
		{
			for (Element aChild : xmlRoot.getChildrenElements("kineticFactor")) {
				iSolute = aSim.getSoluteIndex(aChild.getAttributeValue("solute"));

				// Create and initialise the instance
				_kineticFactor[iFactor] = (IsKineticFactor) (new XMLParser(aChild))
				.instanceCreator("simulator.reaction.kinetic");
				_kineticFactor[iFactor].init(aChild);
				_soluteFactor[iFactor] = iSolute;
				iFactor++;
			}

			_kineticParam = ExtraMath.newDoubleArray(getTotalParam());
			_kineticParam[0] = _muMax;

			paramIndex = 1;
			iFactor = 0;
			for (Element aChild : xmlRoot.getChildrenElements("kineticFactor")) {
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
	 * \brief Use the reaction class to fill the parameters fields of the agent
	 * 
	 * Use the reaction class to fill the parameters fields of the agent. Uses information in the agent and protocol file to achieve this
	 * 
	 * @param anAgent	The ActiveAgent which parameters are being populated
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aReactionRoot	The XML object containing the definition of one reaction in the protocol file
	 * @see Simulator.createReaction()
	 */
	@Override
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot) {
		// Call the init of the parent class (populate yield arrays)
		super.initFromAgent(anAgent, aSim, aReactionRoot);

		anAgent.reactionKinetic[reactionIndex] = ExtraMath.newDoubleArray(getTotalParam());

		// Set muMax
		unit = new StringBuffer("");
		value = aReactionRoot.getParamSuchDbl("kinetic", "muMax", unit);
		double muMax = value*UnitConverter.time(unit.toString());
		anAgent.reactionKinetic[reactionIndex][0] = muMax;

		// Set parameters for each kinetic factor
		paramIndex = 1;
		for (Element aChild : aReactionRoot.getChildrenElements("kineticFactor")) {
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
	 * \brief Update the Marginal Mu data matrix
	 * 
	 * Update the Marginal Mu data matrix
	 * 
	 * @param s	Temporary container for solute concentration 
	 */
	@Override
	public void updateMarginalMu(Double[] s)
	{
		int soluteIndex;
		for (int iFactor = 0; iFactor < _soluteFactor.length; iFactor++)
		{
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1)
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0);
			else
				marginalMu[iFactor] = _kineticFactor[iFactor]
				                                     .kineticValue(s[_soluteFactor[iFactor]]);
		}
	}


	/**
	 * \brief Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * @param s	The concentration of solute locally observed
	 * @param mass	Mass of the catalyst (cell...)
	 * @param tdel	Time
	 */
	@Override
	public void computeUptakeRate(Double[] s, Double mass, Double tdel)
	{
		// First compute specific rate
		computeSpecificGrowthRate(s);
		// Now compute uptake rate
		if(Simulator.isChemostat)
		{
			// TODO Rob 31 Jul 2014: Does this ever get called? Solver_chemostat doesn't seem to use it.
			// (tdel*mass*Dil) is very strange
			for (int iSolute : _mySoluteIndex)
				_uptakeRate[iSolute] = (tdel*mass*Dil) + (mass *_specRate*_soluteYield[iSolute] ) ;
			int iSolute;
			for (int i = 0; i<_soluteFactor.length; i++)
			{
				iSolute = _soluteFactor[i];
				if(iSolute!=-1)	
					_diffUptakeRate[iSolute] =(tdel*mass*Dil) + (mass*marginalDiffMu[i]*_soluteYield[iSolute])  ;	
			}
		}
		else
		{
			for (int iSolute : _mySoluteIndex)
				_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
			int iSolute;
			for (int i = 0; i<_soluteFactor.length; i++)
			{
				iSolute = _soluteFactor[i];
				if(iSolute!=-1)
					_diffUptakeRate[iSolute] = mass*marginalDiffMu[i]*_soluteYield[iSolute];
			}
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
	public Double computeSpecRate(Double[] s)
	{
		Double specRate = _muMax;
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
			specRate *= marginalMu[iFactor];
		return specRate;
	}

	/**
	 * \brief Compute specific growth rate in function to concentrations sent
	 * 
	 * Compute specific growth rate in function to concentrations sent
	 * 
	 * @param s	Array of solute concentration
	 */
	@Override
	public void computeSpecificGrowthRate(Double[] s)
	{
		_specRate = _muMax;
		int soluteIndex;
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
		{
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1)
			{
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0);
				marginalDiffMu[iFactor] = _muMax*_kineticFactor[iFactor].kineticDiff(0.0);
			}
			else
			{
				marginalMu[iFactor] = _kineticFactor[iFactor]
				                                     .kineticValue(s[_soluteFactor[iFactor]]);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]]);
			}
		}

		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++)
		{
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor<_soluteFactor.length; jFactor++)
			{
				if (jFactor!=iFactor)
					marginalDiffMu[jFactor] *= marginalMu[iFactor];
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
	public Double[] getMarginalDiffMu() 
	{
		return marginalDiffMu;
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
	@Override
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
	 * dMUdT = catalyticBiomass*specificGrowthRate*soluteYield.
	 * 
	 * @param S	Temporary container for solute concentration.
	 * @param biomass	Total particle mass in the system which catalyses this reaction.
	 * @return Matrix containing rate of change of each uptake rate with respect to time.
	 */ 
	@Override
	public Matrix calcdMUdT(Matrix S, Double biomass)
	{
		Matrix dMUdT = new Matrix(nSolute, 1, 0.0);
		try
		{
			// No need to initialise this with ExtraMath.newDoubleArray(nI)
			Double[] s = new Double[nSolute];
			for (int i = 0; i < nSolute; i++)
			{
				s[i] = S.get(i,0);
				dMUdT.set(i, 0, _soluteYield[i]);
			}
			updateMarginalMu(s);
			_specRate = computeSpecRate(s);
			dMUdT.timesEquals(_specRate*biomass);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "ReactionFactorWithConstant.calcdMUdT()");
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
	@Override
	public Double[] computeMarginalDiffMu(Double[] s)
	{
		int soluteIndex;
		for (int iFactor = 0; iFactor < _soluteFactor.length; iFactor++)
		{
			soluteIndex = _soluteFactor[iFactor];
			if (soluteIndex==-1)
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0);
			else
			{
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]]);
			}
			for (int jFactor = 0; jFactor < _soluteFactor.length; jFactor++)
			{
				if (jFactor!=iFactor)
					marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
		return marginalDiffMu;
	}
	
	/* __________________ Methods called by the agents ___________________ */
	
	/**
	 * TODO Check that plFitness should be initialised as 0 rather than as 1
	 * 
	 * @param anAgent
	 * @return the marginal growth rate (i.e the specific growth rate times the
	 * mass of the particle which is mediating this reaction)
	 */
	@Override
	public Double computeMassGrowthRate(ActiveAgent anAgent)
	{
		Double plFitness = 0.0;
		for (Double cost : setYield(anAgent))
			plFitness -= cost;
		
		//computing _specRate
		computeSpecificGrowthRate(anAgent);
		
		return _specRate * plFitness * anAgent.getParticleMass(_catalystIndex);
	}

	@Override
	public Double computeSpecGrowthRate(ActiveAgent anAgent)
	{
		Double plFitness = 1.0;
		for (Double cost : setYield(anAgent))
			plFitness -= cost;
		
		//computing _specRate
		computeSpecificGrowthRate(anAgent);
		
		return _specRate*plFitness;
	}
	
	/**
	 * Compute specific growth rate in function to concentrations sent.
	 * 
	 * @param s Double[] array of solute concentrations.
	 * @param anAgent Parameters used are those defined for THIS agent
	 */
	@Override
	public void computeSpecificGrowthRate(Double[] s, ActiveAgent anAgent)
	{
		Double[] kineticParam = anAgent.reactionKinetic[reactionIndex];
		
		// First multiplier is muMax
		_specRate = kineticParam[0];
		paramIndex = 1;
		
		// Compute contribution of each limiting solute
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			if (_soluteFactor[iFactor]==-1) { //meaning, if there's no such solute
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(0.0, kineticParam,
						paramIndex);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(0.0, kineticParam, paramIndex);
				paramIndex += _kineticFactor[iFactor].nParam;
			} else {
				marginalMu[iFactor] = _kineticFactor[iFactor].kineticValue(
						s[_soluteFactor[iFactor]], kineticParam, paramIndex);
				marginalDiffMu[iFactor] = _muMax
				*_kineticFactor[iFactor].kineticDiff(s[_soluteFactor[iFactor]],
						kineticParam, paramIndex);
				paramIndex += _kineticFactor[iFactor].nParam;
			}
		}

		// Finalise the computation
		for (int iFactor = 0; iFactor<_soluteFactor.length; iFactor++) {
			_specRate *= marginalMu[iFactor];
			for (int jFactor = 0; jFactor<_soluteFactor.length; jFactor++) {
				if (jFactor!=iFactor) marginalDiffMu[jFactor] *= marginalMu[iFactor];
			}
		}
	}

	/**
	 * \brief Return the specific reaction rate for a given agent.
	 * 
	 * @param anAgent	Agent to use to determine solute concentration and calculate reaction rate
	 * @see ActiveAgent.grow()
	 * @see Episome.computeRate(EpiBac)
	 */
	@Override
	public void computeSpecificGrowthRate(ActiveAgent anAgent)
	{
		computeSpecificGrowthRate(readConcentrationSeen(anAgent, _soluteList), anAgent);
	}

	/**
	 * Sonia
	 * Here we calculate the fitness cost of each of the plasmid types the host carries. For that we need
	 * the time the plasmid has spent in this cell since its reception or creation (for that we created a specific
	 * field for the plasmid: timeSpentInHost).
	 * @param anAgent
	 */

	public ArrayList<Double> setYield(ActiveAgent anAgent)
	{
		Double initialCost, rateDec, basalCost, plCost, timeSpentInHost;
		int plCopyNum;
		ArrayList<Double> plTotalCosts = new ArrayList<Double>();
		
		if(anAgent instanceof MultiEpiBac)
		{	
			MultiEpiBac anEpiBac = (MultiEpiBac) anAgent;
			//System.out.println("plasmid list size " + anEpiBac._plasmidHosted.size());
			if(anEpiBac.plasmidHosted.size() > 0)
			{	
				for (int pl=0; pl< anEpiBac.plasmidHosted.size(); pl++)
				{
					MultiEpisomeParam plParam = anEpiBac.plasmidHosted.get(pl).getSpeciesParam();
					initialCost = plParam.initialCost;
					//System.out.println("plasmidCost is " + initialCost);
					rateDec = plParam.rateDec;
					basalCost = plParam.basalCost;
					plCopyNum = anEpiBac.plasmidHosted.get(pl).getCopyNumber();
					timeSpentInHost = SimTimer.getCurrentTime()-anEpiBac.plasmidHosted.get(pl).timeSpentInHost;
					//sonia: the cost of a plasmid increases additively as its copy number goes up
					plCost = (initialCost*(Math.exp((-(rateDec*timeSpentInHost))))+ basalCost)*plCopyNum;
					plTotalCosts.add(plCost);					
				}
			}
		}
		return plTotalCosts;	
	}
}
