/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * 
 * 
 */

/**
 * 
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

/**
 * \brief Define the bulk: an object used to define the environment connected to the simulated system
 * 
 * The Bulk is an object used to define the environment connected to the simulated system : this environment can impose concentrations 
 * on the boundaries (constant boundary) or exchange matter through the boundaries (bulk boundary). A  bulk is a perfectly mixed liquid  
 * compartment of the system, usually with a larger size than the extent of the simulated biofilm domain. For example, in a wastewater  
 * reactor the bulk would refer to the liquid volume being treated. The bulk will usually have a fixed volume, but the volume is not 
 * specified here, instead the definition of the computationDomain addresses the bulk compartment volume
 * 
 * @since August 2006
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */
public class Bulk 
{

	/**
	 * Name assigned to this bulk, as specified in the protocol file
	 */
	private String   _name;
	
	/**
	 * The simulation object being used to simulate the conditions specified in the protocol file
	 */
	public Simulator mySim;
	
	/**
	 * Array containing the initial bulk concentration of each solute in the bulk. RJC 22/8/11 - Made public so that Solver_chemostat
	 * can use it to initialise solute concentrations
	 */
	public double[] _bulkValue;
	
	/**
	 * Array of reaction rates for each solute specified in this bulk
	 */
	private double[] _reacRate;
	
	/**
	 * Timestep value assigned to each solute in specified in this bulk
	 */
	Double[]         _dT;
	
	/**
	 * Boolean noting whether the concentration of solutes is constant or varys. Default to constant if not in protocol file
	 */
	private Boolean  _bulkIsConstant = true;
	
	/**
	 * Reactor Dilusion rate of the bulk. Used if isConstant is set to false
	 */
	public  double   _D;
	
	/**
	 * Array of doubles that specify the feed flow concentration of each solute in this simulation (if applicable)
	 */
	public  double[] _sIn;

	/**
	 * Array of booleans that note whether each solute specified in the simulation is present in this bulk
	 */
	private Boolean[] _isInBulk;

	/**
	 * Array of booleans that note whether the concentration of each specified solute varies over time or remains constant
	 */
	public Boolean[] _isConstant;

	/**
	 * Array of sPulse values assigned to each solute, if applicable. Used in pulsing to spike the concentration of a solute at a given rate
	 */
	private double[] _sPulse;
	
	/**
	 * Array of pulseRate values assigned to each solute, if applicable. Sets the rate that the concentration of a solute is spiked if pulsing is employed
	 */
	private double[] _pulseRate;
	
	/**
	 * Array of interval values assigned to each solute, if applicable. Sets the interval at which concentration spiking occurs, if employed
	 */
	private double[] _pulseInterval;
	
	/**
	 * Array that stores the last time in the simulation that particular solute was spiked, if applicable. Used with pulsed concentrations
	 */
	private double[] _lastPulseTime;

	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/

	/**
	 * \brief Creates a bulk compartment object with attributes specified in the protocol file
	 * 
	 * Constructor for the Bulk compartment. For each bulk specified in the protocol file, creates a simulation bulk object with 
	 * the attributes specified in the XML. This includes storing the information for each solute (initial concentrations, reaction rates, 
	 * constant/varying concentration levels, and whether pulsed concentrations are employed) in the retrospective array for access 
	 * in later calculations
	 * 
	 * @param aSim	The simulation object upon which the scenario specified in the protocol file is beinbg run
	 * @param aBulkRoot	The XML tag objects that sub-nodes of the 'Bulk' tag in the protocol file
	 */
	public Bulk(Simulator aSim, XMLParser aBulkRoot) 
	{
		// Solute index is the integer reference to the solute in the simulation solutes dictionary
		int soluteIndex;
		XMLParser parser;

		// Set the simulation object to that initialised in the Simulator class
		mySim = aSim;
		
		// Get the name of this Bulk compartment object
		_name = aBulkRoot.getAttribute("name");
		
		// Check the XML to determine is the concentration of the bulk is assumed to be constant or affected by mass balance
		_bulkIsConstant = aBulkRoot.getParamBool("isConstant");

		if (!_bulkIsConstant) 
		{
			// With bulkIsConstant set to false, the solute concentrations in the bulk vary in time. The protocol file should 
			// specify how the bulk concentration is updated
			String updateType = aBulkRoot.getParam("updateType");
			
			if (updateType != null && updateType.equals("gradient")) 
			{
				//_updateByReaction = false;
				LogFile.writeLog("\t\tUsing gradient-based method for bulk updates.");
			} 
			else
				LogFile.writeLog("\t\tUsing reaction-based method for bulk updates.");
		}

		
		// Array initialisation - store the attributes of each solute in an array
		// Build the list of concentrations for each solute in this bulk
		LinkedList<Element> soluteList = aBulkRoot.buildSetMarkUp("solute");
		_bulkValue = new double[aSim.soluteDic.size()];
		_reacRate = new double[aSim.soluteDic.size()];
		_dT = ExtraMath.newDoubleArray(aSim.soluteDic.size());

		_isConstant = new Boolean[aSim.soluteDic.size()];
		_isInBulk = new Boolean[aSim.soluteDic.size()];
		Arrays.fill(_isInBulk, false);

		_pulseRate = new double[aSim.soluteDic.size()];
		_pulseInterval = new double[aSim.soluteDic.size()];
		_lastPulseTime = new double[aSim.soluteDic.size()];

		_sIn = new double[aSim.soluteDic.size()];
		_sPulse = new double[aSim.soluteDic.size()];
		
		// Parameter D is the Reactor Dilusion rate, used if isConstand is set to false
		_D = aBulkRoot.getParamDbl("D");
	

		// Now iterate through each solute specified in this bulk
		for (Element asoluteMarkUp : soluteList) 
		{
			parser = new XMLParser(asoluteMarkUp);
			// Earlier (in createSimulation method in Simulator) we built a dictionary of solutes in this simulation. Get the index 
			// assigned in that dictionary for this solute name. This is then used as the reference for each of the above arrays
			soluteIndex = aSim.getSoluteIndex(parser.getAttribute("name"));
			
			// Get the Sbulk value for this solute - the initial bulk concentration of this solute
			_bulkValue[soluteIndex] = parser.getParamConc("Sbulk");
			LogFile.writeLog("Setting initial "+_bulkValue[soluteIndex]);
			
			// Note in the boolean array that this solute is present in the bulk if concentration value is over zero
			if (Double.isNaN(_bulkValue[soluteIndex])) _bulkValue[soluteIndex] = 0;
			_isInBulk[soluteIndex]= true;

			// Note in the array whether the concentration of this solute is set to vary over time (as determined earlier in this method
			_isConstant[soluteIndex] = _bulkIsConstant;

			// However each solute can have its own isConstant value - check whether this is the case, and if so determine if the 
			// concentration of this solute is variable. Even if the bulk is set to isConstant, the solute can be set to be varied and 
			// thus override the compartment
			String isconst = parser.getParam("isConstant");
			if (isconst != null)
				_isConstant[soluteIndex] = Boolean.valueOf(isconst);
			
			// Get the concentration of each solute in the feed flow (for varying concentrations)
			_sIn[soluteIndex] = parser.getParamDbl("Sin");
			if (Double.isNaN(_sIn[soluteIndex])) _sIn[soluteIndex] = 0;

			// Now check for pulses. The pulseRate parameters can be used to periodically spike the concentration of a solute to that
			// specified in Spulse parameter at a given time
			
			// for pulses, set pulse concentration and interval
			_sPulse[soluteIndex] = parser.getParamDbl("Spulse");
			if (Double.isNaN(_sPulse[soluteIndex])) _sPulse[soluteIndex] = 0;

			// pulse interval will be infinite if the rate is zero
			_pulseInterval[soluteIndex] = Double.MAX_VALUE;
			_lastPulseTime[soluteIndex] = 0.;
			_pulseRate[soluteIndex] = parser.getParamDbl("pulseRate");
			if (!Double.isNaN(_pulseRate[soluteIndex]) && _pulseRate[soluteIndex]!=0.)
				_pulseInterval[soluteIndex] = 1./_pulseRate[soluteIndex];
			
		}
	}
                  
	/**
	 * \brief Updates the bulk solute concentration if the bulk is not constant
	 * 
	 * Updates the bulk solute concentration if the bulk is not constant.
	 * bvm note 15.12.08: modified routine to allow input of solute grid.
	 * bvm note 13.07.09: removed passing-in of 'implicit' flag, and instead set the bulk update method via the XML protocol file
	 * 
	 * @param soluteGrid	Grid of all solutes in the simulated system
	 * @param reacGrid	Grid of all reactions in the simulated system
	 * @param timeStep	Internal timestep used to update the simulation environment
	 */
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
	//
	/**
	 * \brief Method to update Bulk solute concentrations during a chemostat setup
	 * 
	 * Method to update Bulk solute concentrations during a chemostat setup
	 * 
	 * @author Sonia Martins 190210
	 * 
	 * @param allSol	Grid of all solutes in the simulated system
	 * @param reacGrid	Grid of all reactions in the simulated system
	 */
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
			message += reacGrid[iGrid].gridName+" [";
			message += ExtraMath.toString(_bulkValue[iGrid], true)+" g/L]\n";
		}

		LogFile.writeLog("Bulk update: "+ message);
	}

	/**
	 * \brief Update bulk concentration by reaction, determining reaction rate seen from reaction compartments
	 * 
	 * Update bulk concentration by reaction, determining reaction rate seen from reaction compartments
	 * 
	 * @param reacGrid	An array of uptake-rates grids in g.L-1.h-1
	 * @param timeStep	Internal timestep used to update the simulation environment
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
	 * \brief Update bulk concentration on the basis of the flow passed through the interface with the bulk compartment
	 * 
	 * Update bulk concentration on the basis of the flow passed through the interface with the bulk compartment
	 * 
	 * IMPORTANT: this method should ONLY be used to test the updateBulkByReaction routine with FLAT biofilms because this routine 
	 * does not treat the gradient through the interface correctly.
	 * DO NOT USE FOR GENERAL SIMULATIONS
	 * 
	 * @param soluteGrid	Grid of all solutes in the simulated system
	 * @param timeStep	Internal timestep used to update the simulation environment
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
	 * \brief Compute massic flow
	 * 
	 * Compute massic flow : g.m-2.time-1
	 * 
	 * @param aSG	a solute grid for a particular solute
	 * @param aDC	Position on the grid
	 * @param flow	Flow direction
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

	/**
	 * \brief Determine if a particular solute is in the bulk
	 * 
	 * Determine if a particular solute is in the bulk. Uses the index of this solute in the simulation dictionary
	 * 
	 * @param soluteIndex	Index of solute in the simulation dictionary
	 * @return	Boolean noting whether this solute is in the bulk (true) or not (false)
	 */
	public Boolean contains(int soluteIndex) 
	{
		return _isInBulk[soluteIndex];
	}

	/**
	 * \brief Get the value of a particular solute in the bulk
	 * 
	 * Get the value of a particular solute in the bulk. Uses the index of this solute in the simulation dictionary
	 * 
	 * @param soluteIndex	Index of solute in the simulation dictionary
	 * @return	Level of this particular solute in the bulk
	 */
	public double getValue(int soluteIndex) {
		return _bulkValue[soluteIndex];
	}

	/**
	 * \brief Set the value of a particular solute in the bulk
	 * 
	 * Set the value of a particular solute in the bulk. Uses the index of this solute in the simulation dictionary
	 * 
	 * @param soluteIndex	Index of solute in the simulation dictionary
	 * @param value	Level at which to set the solute level
	 */
	public void setValue(int soluteIndex, double value) {
		_bulkValue[soluteIndex] = value;
	}

	/**
	 * \brief Get the name of this bulk
	 * 
	 * Get the name of this bulk
	 * 
	 * @return	String containing the name of this bulk
	 */
	public String getName() {
		return _name;
	}

	/**
	 * \brief Return the time constraint of the bulk
	 * 
	 * Return the time constraint of the bulk
	 * 
	 * @return	Double containing the time constraint of this bulk
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
	 * \brief Writes a description of the bulk in the result file
	 * 
	 * Writes a description of the bulk in the result file
	 * 
	 * @param buffer	Buffer to which simulation results are being written
	 * @throws Exception	Exception thrown if this buffer cannot be written to
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
