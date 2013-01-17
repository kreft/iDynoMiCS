/**
 * Project iDynoMicS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Keep track of time step and time course
 * 
 */

/**
 * 
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

package idyno;

import simulator.World;

import utils.LogFile;
import utils.XMLParser;

public class SimTimer {

	// Number of iterations
	private static int     _nIter;
	// Size of a time step, in hour
	private static double  _dT, _dTMax, _dTMin;
	private static double[] _oldStep;
	// Store simulation time, in hours
	private static double  _now;
	// Time for the beginning and the end of the simulation, in hours
	private static double  _endOfSimulation;
	private static boolean isAdaptive;

	// _________________ CONSTRUCTOR _____________________________ //

	/**
	 * XML-based constructor
	 * @param localRoot
	 */
	public SimTimer(XMLParser localRoot) {
		XMLParser parser = new XMLParser(localRoot.getChildElement("timeStep"));

		// Set to zero all counters
		reset();

		// Initialize the timer
		_endOfSimulation = parser.getParamTime("endOfSimulation");
		_dT = parser.getParamTime("timeStepIni");

		isAdaptive = parser.getParamBool("adaptive");
		if (isAdaptive) {
			LogFile.writeLog("Using adaptive time stepping.");
			_dTMax = parser.getParamTime("timeStepMax");
			_dTMin = parser.getParamTime("timeStepMin");
			_oldStep = new double[10];
			for(int i=0;i<10;i++)
				_oldStep[i]=_dT;			
		}
	}

	/**
	 * Computes the current simulation time step and increments the simulation
	 * clock.
	 */
	public static void applyTimeStep() {
		_now += _dT;
		_nIter++;
	}

	public static void updateTimeStep(World aWorld) {

		if (!isAdaptive) return;

		double tOpt, newDeltaT;

		tOpt = aWorld.getBulkTimeConstraint();

		if (Double.isInfinite(tOpt)|Double.isNaN(tOpt)) {
			return;
		}

		// constrain value between min and max limits
		newDeltaT = Math.min(Math.max(tOpt, _dTMin), _dTMax);

		if (newDeltaT <= _dT) {
			// if there is a drop in dt, just use the new value
			// (also do this if it is not changing b/c it is at the minimum)
			_dT = newDeltaT;

			// in this case, we also need to re-populate the saved steps
			// so that we don't use too-large values to raise the step again
			for (int i=0;i<10;i++)
				_oldStep[i] = _dT;
		} else {
			// otherwise gradually update the timestep in case of increase

			// make the new option the average of 10 previous steps
			for (int i=1;i<10;i++)
				_oldStep[i]=_oldStep[i-1];
			_oldStep[0] = newDeltaT;
			newDeltaT = utils.ExtraMath.average(_oldStep);

			// again make sure the step isn't too small or too large
			_dT = _oldStep[0] = Math.min(Math.max(newDeltaT,_dTMin),_dTMax);
		}

		// make the step into a nicer number
		_dT = ((int)Math.floor(_dT/(_dTMin/10)))*(_dTMin/10);

		LogFile.writeLog("TimeStep "+_dT+" ("+tOpt+")");
	}

	/**
	 * Sets time step and counters to zero.
	 */
	public static void reset() {
		_now = 0;
		_nIter = 0;
	}

	// ____________________ Get & Set functions ___________________ //
	/**
	 * @return the time spent since the beginning of the simulation (in hour)
	 */
	public static double getCurrentTime() {
		return _now;
	}

	public static int getCurrentIter() {
		return _nIter;
	}

	/**
	 * @return the time spent between two simulation steps (in hour)
	 */
	public static double getCurrentTimeStep() {
		return _dT;
	}

	public static void setCurrentTimeStep(double dt) {
		_dT = dt;
	}

	public void setTimerState(String infoFile) {
		XMLParser fileRoot = new XMLParser(infoFile);
		XMLParser localRoot = new XMLParser(fileRoot.getChildElement("simulation"));

		_nIter = Integer.parseInt(localRoot.getAttribute("iterate"));
		_now = Double.parseDouble(localRoot.getAttribute("time"));

	}

	public static boolean simIsFinished() {
		return (_now>=_endOfSimulation);
	}

	public static boolean isDuringNextStep(double aDate) {
		return aDate>=_now&&aDate<_now+_dT;
	}
}
