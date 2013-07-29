/**
 * \package simulator.detachment
 * 
 * \brief Package of classes that capture detachment of agents from the biomass
 * 
 * Package of classes that capture detachment of agents from the biomass. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.detachment;

import simulator.AgentContainer;
import simulator.agent.LocatedGroup;
import simulator.Simulator;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief SolGrad detachment method. Loop through solutes to find the maximum gradient and use this as a multiplier to calculate detachment rate
 * 
 * The mark-up defines the erosion forces that act on the biofilm surface. Detachment works by removing a layer of biomass based on the 
 * detachment speed and the timestep, with the detachment speed calculated via one of the given forms. This class captures the Proportional 
 * detachment method.
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology, University of Birmingham, UK 
 *
 */
public class DS_SolGrad extends LevelSet 
{
	/**
	 * Constant parameter used to determine the strength of detachment.
	 */
	private double kDet;
	
	/**
	 * Maximum thickness that the biofilm may reach
	 */
	private double maxTh;
	
	private double answer = 0;
	
	private double bulkValue;

	/**
	 * \brief Initialise the object by reading attributes from associated agent grid and XML protocol file
	 * 
	 * Initialise the object by reading attributes from associated agent grid and XML protocol file
	 * 
	 * @param anAgentGrid	Associated grid of agents
	 * @param root	XML tag containing information related to this detachment mechanism
	 */
	public void init(AgentContainer anAgentGrid, XMLParser root) {
		super.init(anAgentGrid, root);
		// kDet has units of: um2.hr-1
		// this gives speed in um.hr-1
		kDet = root.getParamDbl("kDet");
		double value=root.getParamDbl("maxTh");
		maxTh=(Double.isNaN(value)? Double.POSITIVE_INFINITY:value);
	}

	/**
	 *\brief Calculate and return the local detachment speed using this detachment method
	 *
	 * Calculate and return the local detachment speed using this detachment method
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aGroup	Located group for which the local detachment speed is being determined
	 * @return Double stating local detachment speed for this group
	 *
	 */
	protected double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim) 
	{
		// If group is above threshold, definitely remove
		if (aGroup.cc.x>maxTh) {
			aSim.continueRunning = false;
			LogFile.writeLog("Maximum threshold "+maxTh);
			LogFile.writeLog("Simulation halted as threshold crossed at"+aGroup.cc);
			return Double.MAX_VALUE;
		} else 
		{
			// If group below threshold, loop through solutes in order to find the
			// maximum gradient at this point. Note that we don't have to worry about
			// pressure (which is considered a solute in iDynoMiCS), because its bulk
			// value should be 0
			bulkValue = aSim.world.getMaxBulkValue(0);
			if (bulkValue>0)
				if (aSim.is3D){
					answer = aSim.soluteList[0].getGradient(aGroup.cc).norm()/bulkValue;
				}else {
					// in 2D we count the y-direction twice
					answer = aSim.soluteList[0].getGradient2D(aGroup.cc).norm()/bulkValue;
				}

			return kDet*answer;
		}
	}

}