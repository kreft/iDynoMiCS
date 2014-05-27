/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store parameters for these agent types. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.agent.LocatedAgent;

import utils.ExtraMath;

/**
 * \brief Creates an object of the Bacterium EPS Species: a bacterium with a permanent hydrolytic activity of its EPS capsule
 * 
 * The BactEPS species behaves just like the Bacterium species, with one small change: capsular EPS is excreted continuously at a 
 * defined rate rather than during discrete events like for a Bacterium. During excretion, an agent will distribute any excreted EPS 
 * to neighbouring EPS agents of that type, or will create a new EPS agent if none are nearby. This new functionality requires 
 * specification of one additional parameter: the hydrolysis rate kHyd, which controls the rate at which capsular EPS is excreted to 
 * neighboring EPS particles.
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 *
 */
public class BactEPS extends Bacterium 
{
	/**
	 * \brief Constructor used to generate progenitor and initialise an object to store relevant parameters
	 * 
	 * Constructor used to generate progenitor and initialise an object to store relevant parameters
	 */
	public BactEPS() {
		super();
		_speciesParam = new BactEPSParam();
	}

	/**
	 * \brief Called at each time step of the simulation to compute mass growth and update radius, mass, and volume. In this case also manages EPS hydrolysis
	 * 
	 * Called at each time step of the simulation (under the control of the method Step of the class Agent) to compute mass growth 
	 * and update radius, mass, and volume. Also determines whether the agent has reached the size at which it must divide, and 
	 * monitors agent death. In this case also manages EPS hydrolysis
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
	 * \brief Manages EPS hydrolyse by the agent and distribution between neighbouring agents
	 * 
	 * Manages EPS hydrolyse by the agent and distribution between neighbouring agents
	 */
	public void manageEPS() 
	{

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
			double epsMass = particleMass[epsIndex]*(deltaM);
			particleMass[epsIndex] *= 1-deltaM;
			updateSize();
						
			// Distribute to your neighbours
			nEPS = _myNeighbors.size();			
			for (int iNb = 0; iNb<nEPS; iNb++) {
				aNb = _myNeighbors.removeFirst();
				aNb.particleMass[epsIndex] += epsMass/nEPS;
				aNb.updateSize();
			}			
		}
		
		/* Guard against too big bound EPS _______________________________ */
		if (_volume/_totalVolume<(1-getSpeciesParam().epsMax)) {
			double ratio = ExtraMath.getUniRand(.6, .9);
			excreteEPS(ratio);			
		}
	}

	/**
	 * \brief Return the set of parameters that is associated with the object of this species
	 * 
	 * Return the set of parameters that is associated with the object of this species
	 * 
	 * @return Object of BactEPSParam that stores the parameters associated with this species
	 */
	public BactEPSParam getSpeciesParam() {
		return (BactEPSParam) _speciesParam;
	}
}
