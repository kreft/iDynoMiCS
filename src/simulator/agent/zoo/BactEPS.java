/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ___________________________________________________________________________
 * BactEPS : a bacterium with a permanent hydrolytic activity of its EPS capsule 
 * 
 */

/**
 * @since Feb 2008
 * @version 1.0
  * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * 
 */

package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.agent.LocatedAgent;

import utils.ExtraMath;

public class BactEPS extends Bacterium {

	public BactEPS() {
		super();
		_speciesParam = new BactEPSParam();
	}

	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls
	 */
	protected void internalStep() {
		// Compute mass growth over all compartments
		grow();
		updateSize();

		// test if the EPS capsule has to be excreted
		manageEPS();

		// Divide if you have to
		if (willDivide()) divide();

		// Die if you have to
		if (willDie()) die(true);
	}

	/**
	 * Alternative behaviour for EPS slime
	 */
	public void manageEPS() {

		LocatedAgent aNb;
		double deltaM;int nEPS;

		if (!_hasEps) { return; }
		int epsIndex = particleMass.length-1;	
		if (particleMass[epsIndex]==0) return;
		

		
		/* At this line it is sure you have some EPS to hydrolyse __________*/
		
		//Part of the capsule to hydrolyse
		deltaM = 1-Math.exp(-getSpeciesParam().kHyd*SimTimer.getCurrentTimeStep());		
		
		//List all close EPS particles of your EPS species		
		findCloseSiblings(_epsSpecies.speciesIndex);		
				
		if (_myNeighbors.isEmpty()) {
			// Create an EPS agent and send him part of the capsule		
			excreteEPS(deltaM);
			updateSize();
		} else {
			// Update your mass and size			
			double value = particleMass[epsIndex]*(deltaM);
			particleMass[epsIndex] *= 1-deltaM;
			updateSize();
						
			// Distribute to your neighbours
			nEPS = _myNeighbors.size();			
			for (int iNb = 0; iNb<nEPS; iNb++) {
				aNb = _myNeighbors.removeFirst();
				aNb.particleMass[epsIndex] += value/nEPS;
				aNb.updateSize();
			}			
		}
		
		/* Guard against too big bound EPS _______________________________ */
		if (_volume/_totalVolume<(1-getSpeciesParam().epsMax)) {
			double ratio = ExtraMath.getUniRand(.6, .9);
			excreteEPS(ratio);			
		}
	}

	public BactEPSParam getSpeciesParam() {
		return (BactEPSParam) _speciesParam;
	}
}
