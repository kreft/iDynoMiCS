package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.SpatialGrid;

/**
 * Laurent Lardon wrote the prototype AgingBac
 * Modified by Jan Kreft (j.kreft@bham.ac.uk) 2009-06-18
 * Modified by Edd Miles
 * Modified by Rob Clegg (rjc096@bham.ac.uk) Oct 2010 - Jan 2011
 * AgingBac is a bacterium with asymmetric or symmetric division driven by parameter alpha
 * In the published version of the code, ToxicAgingBac was used for the case when the damage
 * was toxic, but this has now been merged into AgingBac together with LinearToxicAgingBac
 * and LinearToxicPintAgingBac. The file ToxicAgingBac is now obsolete but kept for backward
 * compatibility.
 */
public class ToxicAgingBac extends AgingBac
{
	public ToxicAgingBac()
	{
		super();
	}
	
	/**
	* An alternate method of calling fitMassOnGrid which takes toxicity into account.
	* Note that this will make the env_State output files incorrect. 
	*/
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex)
	{
		if ( isDead )
			return;

		Double value = (1-this.age)*particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}
	
	/**
	* An alternate method of calling fitMassOnGrid which takes toxicity into account.
	* Note that this will make the env_State output files incorrect.
	*/
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG)
	{
		if ( isDead )
			return;
		Double value = (1-this.age)*_totalMass/aSpG.getVoxelVolume();
		if (Double.isNaN(value) | Double.isInfinite(value))
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}
	
	/**
	 * 
	 */
	@Override
	public void updateGrowthRates()
	{
		for (int i = 0; i < particleMass.length; i++)
			deltaParticle[i] = 0.0;
		
		Double tStep = SimTimer.getCurrentTimeStep();
		Double pAct =  this.particleMass[0];
		Double pDam = this.particleMass[1];
		Double age = this.age;
		Double Mu = allReactions[reactionActive.get(0)].computeSpecGrowthRate(this);
		Double a = allReactions[reactionActive.get(1)].getKinetic()[0];
		/*
		 * Growth
		 */
		deltaParticle[0] += pAct*Math.expm1(tStep*Mu*(1-age));
		_netGrowthRate = pAct*Mu*(1-age);
		/*
		 * Aging & repair
		 */
		Double ageMass = pAct*Math.expm1(-a*tStep);
		deltaParticle[0] += ageMass;
		deltaParticle[1] -= ageMass;
		repair();
		/*
		 * Removal
		 */
		if (allReactions.length > 3)
		{
			Double k2 = allReactions[3].getKinetic()[0];
			Double k4 = allReactions[4].getKinetic()[0];
			deltaParticle[0] += pAct * Math.expm1( - k2 * tStep );
			deltaParticle[1] += pDam * Math.expm1( - k4 * tStep );
			_netGrowthRate -= pAct * k2;
			_netGrowthRate -= pDam * k4;
		}
		_netVolumeRate = _netGrowthRate/getSpeciesParam().particleDensity[0];
	}
}
	
