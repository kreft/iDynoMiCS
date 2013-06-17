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

import utils.XMLParser;

@Deprecated
/**
 * \brief Modelled First-Order reactions. This class is deprecated
 * 
 * Modelled First-Order reactions. This class is deprecated
 * 
 * @author Andreas D�tsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class FirstOrder extends Reaction {

	// Serial version used for the serialisation of the class
	private static final long   serialVersionUID = 1L;

	private double              _k;
	private static int          iSolute;

	/* _______________________ CONSTRUCTOR _________________________ */
	public void init(Simulator aSim, XMLParser aReactionRoot) {
		super.init(aSim, aReactionRoot);
		_k = aReactionRoot.getParamTime("k");	

		_kineticParam = new double[1];
		_kineticParam[0] = _k;
	}

	/**
	 */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser aReactionRoot) {
		// Call the init of the parent class (populate yield arrays)
		super.initFromAgent(anAgent, aSim, aReactionRoot);


		anAgent.reactionKinetic[reactionIndex] = new double[1];
		anAgent.reactionKinetic[reactionIndex][0] = aReactionRoot.getParamTime("k");

	}

	/* __________________ METHODS _________________________ */

	public void computeUptakeRate(double[] s, ActiveAgent anAgent) {
		_specRate = anAgent.reactionKinetic[reactionIndex][0];
		// Now compute uptake rate and its derivative for each solute
		for (int i = 0; i<_mySoluteIndex.length; i++) {
			iSolute = _mySoluteIndex[i];
			_uptakeRate[iSolute] = anAgent.particleMass[_catalystIndex]*_specRate
			        *anAgent.soluteYield[reactionIndex][iSolute];
			_diffUptakeRate[iSolute] = 0;
		}
	}

	/**
     * @param s	Array of concentration
     * @param mass	Concentration of reactant
     */
	public void computeUptakeRate(double[] s, double mass) {

		_specRate = _k;
		// Now compute uptake rate and its derivative for each solute
		for (int i = 0; i<_mySoluteIndex.length; i++) {
			iSolute = _mySoluteIndex[i];
			_uptakeRate[iSolute] = mass*_specRate*this._soluteYield[iSolute];
			_diffUptakeRate[iSolute] = 0;
		}
	}
	
	/**
     * \brief Return the specific reaction rate
     * 
     * Return the specific reaction rate
     * @param anAgent	The agent
     * @deprecated
     */
	public void computeSpecificGrowthRate(ActiveAgent anAgent) {
		_specRate = anAgent.reactionKinetic[reactionIndex][0];
	}
	
	/**
	 * @param s	Double array
	 * @deprecated
	 * 
	 */
	public void computeSpecificGrowthRate(double[] s) {
		_specRate = this._k;
	}
	
	/**
     * Compute specific growth rate in fonction to concentrations sent
     * @param s	array of solute concentration
     * @param anAgent Parameters used are those defined for THIS agent
     * @deprecated
     */
	public void computeSpecificGrowthRate(double[] s, ActiveAgent anAgent) {
		_specRate = anAgent.reactionKinetic[reactionIndex][0];
	}

	/**
     * return the marginal growth rate (i.e the specific growth rate times the
     * mass of the particle who is mediating this reaction)
     * @param anAgent
     * @return
     */
	public double computeMassGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate*anAgent.getParticleMass(_catalystIndex);
	}

	public double computeSpecGrowthRate(ActiveAgent anAgent) {
		computeSpecificGrowthRate(anAgent);
		return _specRate;
	}

	@Override
	/**
	 * @param s	Double array
	 * @param conc	Concentration factor
	 * @param h	Double
	 * @deprecated
	 * 
	 */
	public void computeUptakeRate(double[] s, double conc, double h) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Matrix calcdMUdS(Matrix S, double biomass) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Matrix calcdMUdT(Matrix S, double biomass) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] computeMarginalDiffMu(double[] s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double computeSpecRate(double[] s) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void updateMarginalMu(double[] s) {
		// TODO Auto-generated method stub
		
	}

	/**
     * Add the contribution of this agent on the reaction grid and the diff
     * reaction grid
     * @see : Reaction.applyReactionCA() and Reaction.applyReactionIbM()
     */

}