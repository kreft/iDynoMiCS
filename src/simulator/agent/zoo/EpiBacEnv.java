package simulator.agent.zoo;

import java.util.*;

import simulator.Simulator;
import simulator.geometry.FluctEnv;;

/**
 * This class describes the different behaviour a specific agent can undergo according to its
 * environment sensibility.
 * @author SoniaMartins
 *
 */

public class EpiBacEnv extends MultiEpiBac{
	
	// Serial version used for the serialisation of the class
	private static final long   serialVersionUID = 1L;

	public boolean died = false ;
    
	/**
	 * Constructor
	 */
	public EpiBacEnv() {
		super();
		_speciesParam = new EpiBacEnvParam();
	}
	
	
public void internalStep() {
		
		
		// Check if some plasmid has a null copy number and remove it if necessary
		checkMissingPlasmid();
		
		// Compute mass growth over all compartments
		grow();
		
		//conjugate();
	
		// test if the EPS capsule has to be excreted
		updateSize();	
		manageEPS();
        
		//assess which environment is present and behave accordingly
		if (Simulator.isFluctEnv){
		StressEffect();
		}
		
		// Divide if you have to
		if (willDivide()) divide();

		// Die if you have to
		if (willDie()) die(true);

	}
		
	
	
	public EpiBacEnvParam getSpeciesParam() {
		return (EpiBacEnvParam) _speciesParam;
	}
	
	/**
	 * The environment itself is set up in another class - FluctEnv.
	 * Here, we deal with the response of the cell to the stress to which it is
	 * sensible to. 
	 * If the cell is sensible to this stress it will die with a certain probability (defined in the protocol file)
	 * given by the field envProbDie.
	 * 
	 */
	public void StressEffect(){
		
		Random num = new Random();
		double rand = num.nextDouble();

		
		ArrayList<String> envSens = this.getSpeciesParam().envSensitivity;
		
		//current environment set by the Deterministic version of FluctEnv
		String currentEnv = FluctEnv.envStatus;
		
		boolean isSensible = false;
	
		//assessing whether the agent is sensible to the present environment
		for (int r=0; r< envSens.size(); r++){
			if(currentEnv.equals(envSens.get(r))){
				isSensible = true;
				break;
			}
		}
	
		//retrieving variable values to kill the agent with a given probability
		HashMap <String, Double> probDieMap = this.getSpeciesParam().envProbDie;
		double probDie=0;
		
		if(isSensible){
				probDie = probDieMap.get(currentEnv);
				if (rand <= probDie){
					this.die(true);
					died=true;
				}else{
				died=false;
			}
		}
		
		//sonia:
		//modify this condition to a more complicated dependence on the environment duration/"intensity"

	}


	/**public String writeOutput() {
		super.writeOutput();	
		//tempString.delete(tempString.length()-2, tempString.length());
				
		tempString.append("diedOfStress,");
		if(died) System.out.println("THIS AGENT DIED DUE TO STRESS EFFECT");
		if(died)tempString.append("1");
		else tempString.append("0");
		tempString.append(";\n");
		
		return tempString.toString();
	}*/
}


