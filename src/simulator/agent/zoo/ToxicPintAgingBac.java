package simulator.agent.zoo;

import utils.ExtraMath;
import idyno.SimTimer;

public class ToxicPintAgingBac extends ToxicAgingBac{

	public ToxicPintAgingBac() {
		super();
	}
	
	@Override
	public boolean willDivide() {

		// this ensures that the checks for when to divide don't occur too often;
		// at most they will occur at the rate of AGENTTIMESTEP
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if (_timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep())
			return false;

		// at this point we will actually check whether to divide
		_timeSinceLastDivisionCheck = 0.0;

		return ExtraMath.radiusOfASphere(particleMass[0]/getSpeciesParam().particleDensity[0]) > this._myDivRadius;
	}
	
}
