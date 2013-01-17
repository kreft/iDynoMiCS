/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * Bulk is an object used to define the environment connected to the simulated 
 * system : this environment can impose concentrations on the boundaries 
 * (constant boundary) or exchange matter through the boundaries (bulk boundary)
 * 
 */

/**
 * @since August 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */

package simulator.geometry;

import idyno.SimTimer;

import java.util.*;

import org.jdom.Element;

import simulator.Simulator;
import simulator.SoluteGrid;
import utils.ExtraMath;

import utils.XMLParser;
import utils.ResultFile;
import utils.LogFile;

public class Bulk {

	private String   _name;
	public Simulator mySim;
	// Rob (22/8/2011): changed _bulkValue to public so that Solver_chemostat can
	// use it to initialise solute concentrations
	//private  double[] _bulkValue;
	public double[] _bulkValue;
	private double[] _reacRate;
	double[]         _dT;
	private Boolean  _bulkIsConstant = true;
	//private Boolean  _updateByReaction = true;

	// Influent load and dilution rate
	//sonia:chemostat
	//changed to public 
	public  double   _D;
	public  double[] _sIn;

	private Boolean[] _isInBulk;

	// to allow setting each solute as constant if needed
	public Boolean[] _isConstant;

	// bvm 30.4.2009
	// variables pertaining to pulsed bulk concentrations
	private double[] _sPulse; // so that pulsed concentration can be different than inflow
	private double[] _pulseRate;
	private double[] _pulseInterval;
	private double[] _lastPulseTime;

	/* __________________ CONSTRUCTOR _________________________________ */

	/**
	 * Constructor based on XML file
	 * @param aBulkRoot
	 */
	public Bulk(Simulator aSim, XMLParser aBulkRoot) {
		int soluteIndex;
		XMLParser parser;

		mySim = aSim;
		_name = aBulkRoot.getAttribute("name");
		_bulkIsConstant = aBulkRoot.getParamBool("isConstant");

		if (!_bulkIsConstant) {
			// set the method used to update the bulk concentration
			String updateType = aBulkRoot.getParam("updateType");
			if (updateType != null && updateType.equals("gradient")) {
				//_updateByReaction = false;
				LogFile.writeLog("\t\tUsing gradient-based method for bulk updates.");
			} else
				LogFile.writeLog("\t\tUsing reaction-based method for bulk updates.");
		}

		// Build the list of concentrations for each solute
		LinkedList<Element> soluteList = aBulkRoot.buildSetMarkUp("solute");
		_bulkValue = new double[aSim.soluteDic.size()];
		_reacRate = new double[aSim.soluteDic.size()];
		_dT = new double[aSim.soluteDic.size()];

		_isConstant = new Boolean[aSim.soluteDic.size()];
		_isInBulk = new Boolean[aSim.soluteDic.size()];
		Arrays.fill(_isInBulk, false);

		_pulseRate = new double[aSim.soluteDic.size()];
		_pulseInterval = new double[aSim.soluteDic.size()];
		_lastPulseTime = new double[aSim.soluteDic.size()];

		// read these values in for both constant and dynamic bulk
		//if (!_bulkIsConstant) {
		_D = aBulkRoot.getParamDbl("D");
		_sIn = new double[aSim.soluteDic.size()];
		_sPulse = new double[aSim.soluteDic.size()];
		//}

		for (Element asoluteMarkUp : soluteList) {
			parser = new XMLParser(asoluteMarkUp);
			soluteIndex = aSim.getSoluteIndex(parser.getAttribute("name"));
			_bulkValue[soluteIndex] = parser.getParamConc("Sbulk");
			LogFile.writeLog("Setting initial "+_bulkValue[soluteIndex]);
			if (Double.isNaN(_bulkValue[soluteIndex])) _bulkValue[soluteIndex] = 0;
			_isInBulk[soluteIndex]= true;

			_isConstant[soluteIndex] = _bulkIsConstant;
			String isconst = parser.getParam("isConstant");
			//if (isconst != null){
			if (isconst != null)
				_isConstant[soluteIndex] = Boolean.valueOf(isconst);
			//}
			// dilution of the solute	
			//	if(!_bulkIsConstant){

			//if(!_isConstant[soluteIndex]){
			_sIn[soluteIndex] = parser.getParamDbl("Sin");
			if (Double.isNaN(_sIn[soluteIndex])) _sIn[soluteIndex] = 0;
			//}


			// for pulses, set pulse concentration and interval
			_sPulse[soluteIndex] = parser.getParamDbl("Spulse");
			if (Double.isNaN(_sPulse[soluteIndex])) _sPulse[soluteIndex] = 0;

			// pulse interval will be infinite if the rate is zero
			_pulseInterval[soluteIndex] = Double.MAX_VALUE;
			_lastPulseTime[soluteIndex] = 0.;
			_pulseRate[soluteIndex] = parser.getParamDbl("pulseRate");
			if (!Double.isNaN(_pulseRate[soluteIndex]) && _pulseRate[soluteIndex]!=0.)
				_pulseInterval[soluteIndex] = 1./_pulseRate[soluteIndex];
			//	}

		}
	}

	// bvm note 15.12.08: modified routine to allow input of solute grid
	// bvm note 13.07.09: removed passing-in of 'implicit' flag, and instead
	//                    set the bulk update method via the XML protocol file
	public void updateBulk(SoluteGrid[] soluteGrid, SoluteGrid[] reacGrid, double timeStep) {

		if (_bulkIsConstant) return;

		// THE GRADIENT METHOD SHOULD ONLY BE USED AS A TEST FOR updateBulkByReaction
		// AND THEN ONLY FOR FLAT BIOFILMS. DO NOT USE FOR GENERAL SIMULATIONS.
		//		// now update by normal dilution process
		//		if (_updateByReaction)
		//			updateBulkByReaction(reacGrid, timeStep);
		//		else
		//			updateBulkByGradient(soluteGrid, timeStep);

		// use only the reaction method in general

		//sonia 19.02.2010 - new method to update the bulk in a chemostat setup
		if(Simulator.isChemostat){
			updateChemostatBulk(soluteGrid, reacGrid);
		}else{

			// bvm 30.4.2009: for pulsed bulk concentrations
			for (int i=0; i<_bulkValue.length; i++)
				if (SimTimer.getCurrentTime()-_lastPulseTime[i] >= _pulseInterval[i]) {
					// set bulk values to pulsed concentration
					_bulkValue[i] = _sPulse[i];
					_lastPulseTime[i] = SimTimer.getCurrentTime();
					LogFile.writeLog("Pulsed bulk concentration for "+
							mySim.soluteList[i].getName()+": "+_sPulse[i]);
				}
			updateBulkByReaction(reacGrid, timeStep);
		}
	}


	//sonia 19.02.2010
	//method to update Bulk solute concentrations during a chemostat setup
	public void updateChemostatBulk(SoluteGrid[] allSol, SoluteGrid[] reacGrid){
		String message = "Bulk dynamics \n";	

		for (int iGrid = 0; iGrid<allSol.length; iGrid++) {

			//_reacRate[iGrid] = reacGrid[iGrid].getAverageChemo();
			//_bulkValue[iGrid]=allSol[iGrid].getAverageChemo() ;
			_reacRate[iGrid] = reacGrid[iGrid].grid[0][0][0];
			_bulkValue[iGrid]=allSol[iGrid].grid[0][0][0] ;
		}

		// bvm 30.4.2009: for pulsed bulk concentrations
		for (int i=0; i<_bulkValue.length; i++)
			if (SimTimer.getCurrentTime()-_lastPulseTime[i] >= _pulseInterval[i]) {
				// set bulk values to pulsed concentration
				_bulkValue[i] = _sPulse[i];
				_lastPulseTime[i] = SimTimer.getCurrentTime();
				LogFile.writeLog("Pulsed bulk concentration for "+
						mySim.soluteList[i].getName()+": "+_sPulse[i]);
			}
		for (int iGrid = 0; iGrid<allSol.length; iGrid++) {
			//System.out.println ("BULK SOLUTE VALUES......" + _bulkValue[iGrid]);
			message += reacGrid[iGrid].gridName+" [";
			message += ExtraMath.toString(_bulkValue[iGrid], true)+" g/L]\n";
		}

		LogFile.writeLog("Bulk update: "+ message);
	}

	/**
	 * Update bulk concentration
	 * @param reacGrid : an array of uptake-rates grids in g.L-1.h-1
	 * @param timeStep
	 */
	public void updateBulkByReaction(SoluteGrid[] reacGrid, double timeStep) {

		Domain aDomain;
		double factor, dSdT, oldValue;

		String message = "Bulk dynamics \n";

		// Determine reaction rate seen from reactor compartment
		for (int iGrid = 0; iGrid<reacGrid.length; iGrid++) {
			// we don't do the bulk update for certain solutes
			if (reacGrid[iGrid]==null) continue;
			//if (reacGrid[iGrid].gridName.contentEquals("o2d")) continue;
			if (reacGrid[iGrid].gridName.contentEquals("pressure")) continue;
			//Brian 14.04.2010
			if (_isConstant[iGrid]) continue;

			oldValue = _bulkValue[iGrid];

			aDomain = reacGrid[iGrid].getDomain();

			// Scale-up factor (V_domain to V_reactor scaling in m3/m3)
			factor = (aDomain.length_X*1e-6) * aDomain.specificArea;

			// Compute average reaction rate in g.L-1.h-1
			_reacRate[iGrid] = reacGrid[iGrid].getAverage()*factor;

			// compute the total rate change by dilution and reaction
			dSdT = _D*(_sIn[iGrid]-_bulkValue[iGrid]) + _reacRate[iGrid];

			// now do the actual update of the concentration

			// forward Euler method
			//			_bulkValue[iGrid] = _bulkValue[iGrid]+dSdT*timeStep;

			// backward Euler method
			_bulkValue[iGrid] = (_bulkValue[iGrid] +
					timeStep*(_D*_sIn[iGrid]+_reacRate[iGrid]))
					/(1+timeStep*_D);

			// safety catch on concentration
			if (_bulkValue[iGrid] < 0.) _bulkValue[iGrid] = 0.;

			// finally update the timestep value
			// safety catch if nothing is happening (to avoid infinite dt)
			if (dSdT == 0.)
				_dT[iGrid] = Double.MAX_VALUE;
			else {
				// these times are not used to update, but stored for changing the timestep
				double t1, t2;
				// t1 is the time needed to change the bulk by 5%
				t1 = 0.05*_bulkValue[iGrid]/Math.abs(dSdT);
				// t2 is 5% of the time to return to the influent value
				t2 = 0.05*(_bulkValue[iGrid]-_sIn[iGrid])/dSdT;
				if (t2<0)
					// t2 is only used for the situation that S is returning toward Sin
					_dT[iGrid] = Math.min(t1, -t2);
				else
					_dT[iGrid]=t1;
			}

			// some output
			message += reacGrid[iGrid].gridName+" [";
			message += ExtraMath.toString(oldValue, true)+" -> "
			+ExtraMath.toString(_bulkValue[iGrid], true);
			message += " ("+ExtraMath.toString(dSdT, true)+")";
			message += " step "+_dT[iGrid]+" ]\n";
		}

		LogFile.writeLog("Bulk update: "+message);
	}

	/**
	 * Update bulk concentration on the basis of the flow passed through the
	 * interface with the bulk compartment
	 * 
	 * IMPORTANT: this method should ONLY be used to test the updateBulkByReaction
	 * routine with FLAT biofilms because this routine does not treat the
	 * gradient through the interface correctly.
	 * 
	 * DO NOT USE FOR GENERAL SIMULATIONS
	 * 
	 * @param soluteGrid
	 */
	public void updateBulkByGradient(SoluteGrid[] soluteGrid, double timeStep) {

		double volRate, oldValue, dSdT;
		Domain aDomain;
		LinkedList<DiscreteVector> border;
		ContinuousVector flow;

		String message = "Bulk dynamics \n";

		// Now compute massic flow through the boundary
		for (int iGrid = 0; iGrid<soluteGrid.length; iGrid++) {
			// we don't do the bulk update for certain solutes
			if (soluteGrid[iGrid]==null) continue;
			//if (soluteGrid[iGrid].gridName.contentEquals("o2d")) continue;
			if (soluteGrid[iGrid].gridName.contentEquals("pressure")) continue;
			if (_isConstant[iGrid]) continue;

			oldValue = _bulkValue[iGrid];

			aDomain = (Domain) soluteGrid[iGrid].getDomain();

			// Ask for the interface between bulk and D-R field
			border = aDomain.getBorder();

			// Sum the flow; this yields units of [g.m-2.h-1]
			flow = new ContinuousVector(0, 0, 0);
			for (DiscreteVector aDC : border) {
				computeFlow(soluteGrid[iGrid], aDC, flow);
			}
			// flow.x is a SUMMATION over all points on the border,
			// and so we divide to get the average flux per unit area of the border
			// [g.m-2.h-1]
			// (note that because we only use flow.x, it's as if we have a flat
			//  interface all the time (normal is only in x direction))
			double flux = flow.x/border.size();

			// fraction by which the flux border is larger/smaller than the
			// domain carrier area (need the '1.' to keep it a double);
			// the flux to/from the bulk will be per carrier area, not border area
			double coverage = 1.*border.size()/aDomain._nJ/aDomain._nK;

			// calculate concentration changes from flux and geometry:
			// need to change the flow into volume variation so that
			// it will change the bulk concentration correctly; use specificArea
			// to convert to volume rate and 1e-3 to yield m3 -> L volume units
			// [g.L-1.h-1]
			volRate = flux*coverage*aDomain.specificArea*1e-3;

			// compute the total rate change by dilution and reaction
			dSdT = _D*(_sIn[iGrid]-_bulkValue[iGrid]) + volRate;

			LogFile.writeLog("Value Sin: "+_sIn[iGrid]);
			LogFile.writeLog("Value bulkVal: "+_bulkValue[iGrid]);
			LogFile.writeLog("Value volrate: "+volRate);


			// bvm 12.12.08: set this to use Backward Euler like other routine
			_bulkValue[iGrid] = (_bulkValue[iGrid] +
					timeStep*(_D*_sIn[iGrid]+volRate))
					/(1+timeStep*_D);

			// safety catch on concentration
			if (_bulkValue[iGrid] < 0.) _bulkValue[iGrid] = 0.;

			// finally update the timestep value
			// safety catch if nothing is happening (to avoid infinite dt)
			if (dSdT == 0.)
				_dT[iGrid] = Double.MAX_VALUE;
			else {
				// these times are not used to update, but stored for changing the timestep
				double t1, t2;
				// t1 is the time needed to change the bulk by 5%
				t1 = 0.05*_bulkValue[iGrid]/Math.abs(dSdT);
				// t2 is 5% of the time to return to the influent value
				t2 = 0.05*(_bulkValue[iGrid]-_sIn[iGrid])/dSdT;
				if (t2<0)
					// t2 is only used for the situation that S is returning toward Sin
					_dT[iGrid] = Math.min(t1, -t2);
				else
					_dT[iGrid]=t1;
			}

			// some output
			message += soluteGrid[iGrid].gridName+" [";
			message += ExtraMath.toString(oldValue, true)+" -> "
			+ExtraMath.toString(_bulkValue[iGrid], true);
			message += " ("+ExtraMath.toString(dSdT, true)+")";
			message += " step "+_dT[iGrid]+" ]\n";
		}

		LogFile.writeLog("Bulk update: "+message);
	}

	/**
	 * Compute massic flow : g.m-2.time-1
	 * @param aSG
	 * @param aDC
	 * @param flow
	 */
	public void computeFlow(SoluteGrid aSG, DiscreteVector aDC, ContinuousVector flow) {
		int _i, _j, _k;
		_i = aDC.i;
		_j = aDC.j;
		_k = aDC.k;
		double D = aSG.getDiffusivity(); // this is already in units um2/hour
		double[][][] u = aSG.grid;       // units of fg/um3
		double r = aSG.getResolution();  // units of um

		// factor converts units: from fg/(um2.h) to g/(m2.h)
		double factor = 1e-3;

		// this should really be the product of grad(u) and the interface normal,
		// but that normal doesn't seem very easy to compute
		flow.x += -D*(u[_i+1][_j][_k]-u[_i-1][_j][_k])/(2*r)*factor;
		flow.y += -D*(u[_i][_j+1][_k]-u[_i][_j-1][_k])/(2*r)*factor;
		flow.z += -D*(u[_i][_j][_k+1]-u[_i][_j][_k-1])/(2*r)*factor;
	}

	/* _________________________ GET & SET ________________________________ */

	public Boolean contains(int soluteIndex) {
		return _isInBulk[soluteIndex];
	}

	public double getValue(int soluteIndex) {
		return _bulkValue[soluteIndex];
	}

	public void setValue(int soluteIndex, double value) {
		_bulkValue[soluteIndex] = value;
	}

	public String getName() {
		return _name;
	}

	/**
	 * Send the time constraint of that bulk
	 * @return
	 */
	public double getTimeConstraint() {
		double out = ExtraMath.max(_dT);
		for (int iGrid = 0; iGrid<_dT.length; iGrid++) {
			if (_dT[iGrid]==0) continue;
			out = Math.min(out, Math.abs(_dT[iGrid]));
		}
		if (out==0) out = Double.POSITIVE_INFINITY;

		return out;
	}

	/**
	 * Insert the description of the bulk in the result file
	 * @param buffer
	 */
	public void writeReport(ResultFile buffer) throws Exception {
		StringBuffer text = new StringBuffer();
		String soluteName;

		text.append("<bulk name=\"").append(_name).append("\">\n");

		// Concentration
		// bvm 23.01.09: added units to the output (ARE THEY CORRECT?)
		for (int i = 0; i<_bulkValue.length; i++) {
			soluteName = mySim.soluteDic.get(i);
			text.append("<solute name=\"")
			.append(soluteName)
			.append("\" unit=\"g.L-1\">");
			text.append(_bulkValue[i]).append("</solute>\n");
		}
		buffer.write(text.toString());

		// bvm 23.01.09: added units to the output (ARE THEY CORRECT?)
		text = new StringBuffer();
		for (int i = 0; i<_bulkValue.length; i++) {
			soluteName = mySim.soluteDic.get(i);
			// bvm 04.12.08: modified the XML output to make it valid
			text.append("<uptake_rate name=\"")
			.append(soluteName)
			.append("\" unit=\"g.L-1.hour-1\">");
			text.append(this._reacRate[i]).append("</uptake_rate>\n");
		}
		buffer.write(text.toString());
		buffer.write("</bulk>\n");
	}
}
