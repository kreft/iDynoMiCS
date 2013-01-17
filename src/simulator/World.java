/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * ______________________________________________________
 */

package simulator;

import java.util.*;

import utils.ExtraMath;
import utils.XMLParser;
import simulator.geometry.*;

public class World {

	public ArrayList<Domain> domainList = new ArrayList<Domain>();
	public LinkedList<Bulk>  bulkList   = new LinkedList<Bulk>();
	public double[]          bulkTime   = new double[2];

	/* _____________________ CONSTRUCTOR ________________________ */

	/**
	 * World constructor Defines domain(s) and bulk(s)
	 */

	public void init(Simulator aSim, XMLParser worldRoot) {
		// Create & register the defined bulks
		for (XMLParser aBulkRoot : worldRoot.buildSetParser("bulk"))
			bulkList.add(new Bulk(aSim, aBulkRoot));

		// Create & register the defined domains
		for (XMLParser aDomainRoot : worldRoot.buildSetParser("computationDomain"))
			domainList.add(new Domain(aSim, aDomainRoot));

	}

	/* _________________________ TOOLS _________________________________ */

	public Domain getDomain(String cDName) {
		for (int i = 0; i<domainList.size(); i++) {
			if (domainList.get(i).getName().equals(cDName)) { return (Domain) domainList.get(i); }
		}
		return null;
	}

	//
	public Bulk getBulk(String bulkName) {
		for (Bulk aBulk : bulkList) {
			if (aBulk.getName().equals(bulkName)) { return aBulk; }
		}
		return bulkList.getFirst();
	}

	public boolean containsBulk(String bulkName) {
		for (Bulk aBulk : bulkList) {
			if (aBulk.getName().equals(bulkName)) { return true; }
		}
		return false;
	}

	/**
	 * @return the time needed to change of 100% a bulk concentration
	 */
	public double getBulkTimeConstraint() {
		for(int i=1;i<bulkTime.length;i++)
			bulkTime[i]=bulkTime[i-1];

		bulkTime[0] = bulkList.getFirst().getTimeConstraint();
		for (Bulk aBulk : bulkList)
			bulkTime[0] = Math.min(bulkTime[0], aBulk.getTimeConstraint());

		return bulkTime[0];

	}

	public double[] getAllBulkValue(int soluteIndex) {
		double value[] = new double[bulkList.size()];
		for (int i = 0; i<bulkList.size(); i++) {
			if (bulkList.get(i).contains(soluteIndex))
				value[i] = bulkList.get(i).getValue(soluteIndex);
			else
				value[i] = 0.;
		}
		return value;
	}

	public double getMaxBulkValue(int soluteIndex) {
		// Rob 4/3/11: simplified
		return ExtraMath.max(getAllBulkValue(soluteIndex));
		//double maxValue = bulkList.get(0).getValue(soluteIndex);
		//for (int i = 1; i<bulkList.size(); i++) {
		//	maxValue = Math.max(maxValue, bulkList.get(i).getValue(soluteIndex));
		//}
		//return maxValue;
	}

}
