/**
 * \package reaction
 * \brief Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS
 * 
 * Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.reaction;

import Jama.Matrix;
import simulator.Simulator;
import simulator.agent.*;
import utils.UnitConverter;
import utils.ExtraMath;
import utils.XMLParser;


@Deprecated
/**
 * \brief Deprecated. Was used  to create pathways described by a simple Monod kinetic y1 S1 -> y2 S2+X where µ=muMax*S1/(Ks+S1) 
 * 
 * Deprecated. Was used  to create pathways described by a simple Monod kinetic y1 S1 -> y2 S2+X where µ=muMax*S1/(Ks+S1) the substrate 
 * used by a Monod kinetic is the first one of the solutelist of the current object
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class SingleMonod extends Reaction {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	private Double     _Ks, _muMax;
	private static int iSolute;

	/* ________________________ CONSTRUCTORS ________________________________ */
	public SingleMonod() {
	}

	public void init(Simulator aSim, XMLParser aReactionRoot) {
		// Call the init of the parent class (populate yield arrays)
		super.init(aSim, aReactionRoot);

		StringBuffer unit = new StringBuffer("");
		Double value = aReactionRoot.getParamSuchDbl("kinetic", "muMax",unit);
		_muMax = value*UnitConverter.time(unit.toString());
		
		_Ks = aReactionRoot.getParamSuchDbl("kinetic", "Ks");

		_kineticParam = ExtraMath.newDoubleArray(2);
		_kineticParam[0] = _muMax;
		_kineticParam[1] = _Ks;
	}

	/**
     * Used to initialise reaction parameters of the agent
     */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot) {
		// Call the init of the parent class (populate yield arrays)
		super.initFromAgent(anAgent, aSim, aReactionRoot);

		anAgent.reactionKinetic[reactionIndex] = ExtraMath.newDoubleArray(2);
		
		StringBuffer unit = new StringBuffer("");
		Double value = aReactionRoot.getParamSuchDbl("kinetic", "muMax",unit);		
		anAgent.reactionKinetic[reactionIndex][0] = value*UnitConverter.time(unit.toString());
		
		anAgent.reactionKinetic[reactionIndex][1] = aReactionRoot.getParamSuchDbl("kinetic", "Ks",unit);
	}

	/* ___________________________ METHODS __________________________________ */

	/**
     * Return the specific reaction rate
     * @see ActiveAgent.grow()
     * @see Episome.computeRate(EpiBac)
     */
	public void computeSpecificGrowthRate(ActiveAgent anAgent)
	{
		int localIndex = _mySoluteIndex[0];
		Double[] s = readConcentrationSeen(anAgent, _soluteList);
		_specRate = kineticValue(s[localIndex], anAgent.reactionKinetic[reactionIndex], 1);
	}

	/**
     * \brief Return the specific reaction rate
     * 
     * Return the specific reaction rate
     * @param s	array of solute concentration
     * @param anAgent	The agent
     * @deprecated
     */
	public void computeSpecificGrowthRate(Double[] s, ActiveAgent anAgent) {
		int localIndex = _mySoluteIndex[0];		
		_specRate = kineticValue(s[localIndex], anAgent.reactionKinetic[reactionIndex], 1);			
	}
	
	/**
     * \brief Compute specific growth rate in function of concentrations sent. Parameters used are those defined for the reaction
     * 
     * Compute specific growth rate in function of concentrations sent. Parameters used are those defined for the reaction
     * 
     * @param s	array of solute concentration
     * @deprecated
     */
	public void computeSpecificGrowthRate(Double[] s) {
		int localIndex = _mySoluteIndex[0];
		_specRate = _muMax*s[localIndex]/(_Ks+s[localIndex]);
	}
	
	
	/**
     * return the marginal growth rate (i.e the specific growth rate times the
     * mass of the particle who is mediating this reaction)
     * @param anAgent
     * @return
     */
	public Double computeMassGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate*anAgent.getParticleMass(_catalystIndex);
	}

	public Double computeMassGrowthRate(ActiveAgent anAgent, Double[] reactionKinetic) {
		int localIndex = _mySoluteIndex[0];
		Double[] s = readConcentrationSeen(anAgent, _soluteList);
		_specRate = reactionKinetic[0]*kineticValue(s[localIndex], reactionKinetic, 1);

		return _specRate*anAgent.getParticleMass(_catalystIndex);
	}

	public Double computeSpecGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate;
	}

	/**
     * Add the contribution of this agent on the reaction grid and the diff
     * reaction grid
     * @see : Reaction.applyReactionCA() and Reaction.applyReactionIbM()
     */
	public void computeAllUptakeRate(Double[] s, Double mass, Double[] soluteYield,
	        Double[] kineticParam) {
		int localIndex = _mySoluteIndex[0];
		Double s1 = s[localIndex];

		_specRate = kineticParam[0]*kineticValue(s[localIndex], kineticParam, 1);

		// Now compute uptake rate and its derivative for each solute
		for (int i = 0; i<_mySoluteIndex.length; i++) {
			iSolute = _mySoluteIndex[i];
			_uptakeRate[iSolute] = mass*_specRate*soluteYield[iSolute];
			_diffUptakeRate[iSolute] = 0.0;
		}
		_diffUptakeRate[localIndex] = mass*soluteYield[localIndex]*soluteYield[0]
		        *kineticDiff(s1, kineticParam, 1);
	}

	/**
	 * @param s	Double array
	 * @param mass	Concentration factor
	 * @param tdel	Double
	 * @deprecated
	 * 
	 */
	public void computeUptakeRate(Double[] s, Double mass, Double tdel) {
		int localIndex = _mySoluteIndex[0];
		Double s1 = s[localIndex];

		computeSpecificGrowthRate(s);
		//sonia:chemostat 27.11.09
		if(Simulator.isChemostat){
			
			for (int i = 0; i<_mySoluteIndex.length; i++) {
				iSolute = _mySoluteIndex[i];
				_uptakeRate[iSolute] = (tdel* mass*Dil) + (mass *_specRate*_soluteYield[iSolute] ) ;
				_diffUptakeRate[iSolute] = 0.0;
		}
	
			_diffUptakeRate[localIndex] = (tdel*mass*Dil) + mass*_soluteYield[localIndex]*_soluteYield[0]
			                             *kineticDiff(s1, _kineticParam, 1);
						
		}else{
		// Now compute uptake rate and its derivative for each solute
		for (int i = 0; i<_mySoluteIndex.length; i++) {
			iSolute = _mySoluteIndex[i];
			_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
			_diffUptakeRate[iSolute] = 0.0;
		}
		_diffUptakeRate[localIndex] = mass*_soluteYield[localIndex]*_soluteYield[0]
		        *kineticDiff(s1, _kineticParam, 1);
		}
	}

	public void computeUptakeRate(Double[] s, ActiveAgent anAgent) {
		int localIndex = _mySoluteIndex[0];
		Double s1 = s[localIndex];
		
		// First compute specific rate
		computeSpecificGrowthRate(s, anAgent);

		Double mass = anAgent.particleMass[_catalystIndex];

		// Now compute uptake rate and its derivative for each solute
		for (int i = 0; i<_mySoluteIndex.length; i++) {
			iSolute = _mySoluteIndex[i];
			_uptakeRate[iSolute] = mass*_specRate*_soluteYield[iSolute];
			_diffUptakeRate[iSolute] = 0.0;
		}
		_diffUptakeRate[localIndex] = mass*_soluteYield[localIndex]*_soluteYield[0]
		        *kineticDiff(s1, _kineticParam, 1);
	}
	
	
	
	public Double kineticValue(Double solute, Double[] paramTable, int index)
	{
		return solute/(paramTable[index]+solute);
	}

	public Double kineticValue(Double solute) {
		return solute/(_Ks+solute);
	}

	public Double kineticDiff(Double solute, Double[] paramTable, int index) {
		return paramTable[index]/(ExtraMath.sq(paramTable[index]+solute));
	}

	public Double kineticDiff(Double solute) {
		return _Ks/(ExtraMath.sq(_Ks+solute));
	}

	public Double kineticMax() {
		return 1.0;
	}

	public Double kineticMax(Double[] paramTable, int index) {
		return 1.0;
	}


	@Override
	public Matrix calcdMUdS(Matrix S, Double biomass) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Matrix calcdMUdT(Matrix S, Double biomass)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Double[] computeMarginalDiffMu(Double[] s)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Double computeSpecRate(Double[] s)
	{
		// TODO Auto-generated method stub
		return 0.0;
	}

	@Override
	public void updateMarginalMu(Double[] s)
	{
		// TODO Auto-generated method stub
		
	}



}
