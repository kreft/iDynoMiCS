/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 * 
 */

/**
 * @since Mar 2011
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology, University of Birmingham, UK 
 */

package simulator.detachment;

import simulator.AgentContainer;
import simulator.agent.LocatedGroup;
import simulator.Simulator;
import utils.LogFile;
import utils.XMLParser;

public class DS_SolGrad extends LevelSet {

	private double kDet;
	private double maxTh;
	private double answer = 0;
	private double bulkValue;

	public void init(AgentContainer anAgentGrid, XMLParser root) {
		super.init(anAgentGrid, root);
		// kDet has units of: um2.hr-1
		// this gives speed in um.hr-1
		kDet = root.getParamDbl("kDet");
		double value=root.getParamDbl("maxTh");
		maxTh=(Double.isNaN(value)? Double.POSITIVE_INFINITY:value);
	}

	protected double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim) {
		// If group is above threshold, definitely remove
		if (aGroup.cc.x>maxTh) {
			aSim.continueRunning = false;
			LogFile.writeLog("Maximum threshold "+maxTh);
			LogFile.writeLog("Simulation halted as threshold crossed at"+aGroup.cc);
			return Double.MAX_VALUE;
		} else {
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