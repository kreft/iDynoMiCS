package simulator.agent.zoo;

import idyno.SimTimer;

public class LinearToxicAgingBac extends AgingBac{

	public LinearToxicAgingBac() {
		super();
	}
	
	@Override
	public void updateGrowthRates()
	{
		for (int i = 0; i < particleMass.length; i++)
			deltaParticle[i] = 0.0;
		
		Double tStep = SimTimer.getCurrentTimeStep();
		Double pAct = this.particleMass[0];
		Double pDam = this.particleMass[1];
		Double age = this.age;
		Double Mu = allReactions[reactionActive.get(0)].computeSpecGrowthRate(this);
		Double a = allReactions[2].getKinetic()[0];
		/*
		 * Growth
		 */
		deltaParticle[0] += (1-age)*Mu*tStep;
		_netGrowthRate = (1-age)*Mu;
		/*
		 * Aging
		 */
		Double ageMass = pAct*Math.expm1(-a*tStep);
		deltaParticle[0] += ageMass;
		deltaParticle[1] -= ageMass;
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
