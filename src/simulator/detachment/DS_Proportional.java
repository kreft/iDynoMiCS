package simulator.detachment;

import simulator.AgentContainer;
import simulator.Simulator;
import simulator.agent.LocatedGroup;
import utils.XMLParser;

public class DS_Proportional extends LevelSet {

	private double kDet;
	private double maxTh;

	public void init(AgentContainer anAgentGrid, XMLParser root) {
		super.init(anAgentGrid, root);
		// kDet has units of: fg.um-4.hr-1
		// this gives speed in um.hr-1
		kDet = root.getParamDbl("kDet");
		double value=root.getParamDbl("maxTh");
		maxTh=(Double.isNaN(value)? Double.POSITIVE_INFINITY:value);
	}

	protected double getLocalDetachmentSpeed(LocatedGroup aGroup, Simulator aSim) {
		if (aGroup.cc.x>maxTh) return Double.MAX_VALUE;
		return kDet*aGroup.cc.x;
	}

}