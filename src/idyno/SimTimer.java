/**
 * \package idyno
 * \brief Package of classes used to launch iDynomics
 * 
 * Package of classes used to launch and iDynoMiCS simulation, and to update the package to the latest stable release. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the 
 * following URL  "http://www.cecill.info".
 */
package idyno;

import simulator.World;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Class to create and keep track of the timestep and simulation time course
 * 
 * Class to create and keep track of the timestep and simulation time course
 *
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class SimTimer {

	/**
	 * Number of simulation iterations performed
	 */
	private static int     _nIter;
	
	/**
	 * Time step of the simulation, in hours. Read in from the timeStepIni parameter of the XML protocol file
	 */
	private static double  _dT;
	
	/**
	 * As _dT changes throughout the simulation, statedTimeStep holds the value that was originally read into _dT
	 */
	public double _statedTimeStep;
	
	/**
	 * Boolean flag that notes whether adaptive time steps are being employed. Read in from adaptive parameter in XML protocol file
	 */
	private static boolean isAdaptive;
	
	/**
	 * Maximum time step value (in hours) that can be used by the simulator, when adaptive time steps are employed. Read in from 
	 * timeStepMax parameter in XML protocol file
	 */
	private static double	_dTMax;
	
	/**
	 * Minimum time step value (in hours) that can be used by the simulator, when adaptive time steps are employed.  Read in from 
	 * timeStepMin parameter in XML protocol file
	 */
	private static double	_dTMin;
	
	
	private static Double[] _oldStep;
	
	/**
	 * Simulation time - in hours
	 */
	private static double  _now;
	
	/**
	 * Amount of time for which the simulation should run. Specified in hours
	 */
	private static double  _endOfSimulation;
	
	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/
	
	
	
	/**
	 * \brief Parses simulation time step information from the XML protocol file
	 * 
	 * The timeStep markup within the protocol file specified all parameters for the time-stepping of the simulation. This action 
	 * is controlled by an object of this SimTimer class. There is a choice of using an adaptive or set timestep, specified by setting 
	 * the adaptive parameter to true or false respectively. If an adaptive timestep is used, three timeStep parameters are set to 
	 * control the initial, minimum, and maximum values the timestep may take. The simulation then determines during runtime the timestep 
	 * that should be used. On the other hand, if this is set to false, the simulation runs at one timestep - that set in the parameter 
	 * timeStepIni. A further parameter is read in, endOfSimulation, that specifies when the simulation should end
	 * 
	 * @param localRoot	The SIMULATOR mark-up within the XML file - containing the timestep parameters.
	 */
	public SimTimer(XMLParser localRoot) 
	{
		XMLParser parser = new XMLParser(localRoot.getChildElement("timeStep"));

		// Set to zero all counters
		reset();

		// Initialize the timer
		_endOfSimulation = parser.getParamTime("endOfSimulation");
		_dT = parser.getParamTime("timeStepIni");
		_statedTimeStep = parser.getParamTime("timeStepIni");

		// Now determine if adaptive timesteps are being used
		isAdaptive = parser.getParamBool("adaptive");
		
		if (isAdaptive) 
		{
			LogFile.writeLog("Using adaptive time stepping.");
			_dTMax = parser.getParamTime("timeStepMax");
			_dTMin = parser.getParamTime("timeStepMin");
			_oldStep = ExtraMath.newDoubleArray(10);
		}
	}

	/**
	 * \brief Computes the current simulation time step and increments the simulation clock.
	 * 
	 * Computes the current simulation time step and increments the simulation clock. First method called by the simulation step method
	 */
	public static void applyTimeStep() 
	{
		_now += _dT;
		_nIter++;
	}

	/**
	 * \brief If adaptive timesteps are enabled, updates the timestep for a given simulation world
	 * 
	 * If adaptive timesteps are enabled, updates the timestep for a given simulation world
	 * 
	 * @param aWorld	The simulation world being modelled under an adaptive timestep
	 */
	public static void updateTimeStep(World aWorld)
	{
		if (!isAdaptive) return;
		
		double tOpt, newDeltaT;
		
		tOpt = aWorld.getBulkTimeConstraint();
		
		if (Double.isInfinite(tOpt)|Double.isNaN(tOpt))
			return;
		
		// constrain value between min and max limits
		newDeltaT = Math.min(Math.max(tOpt, _dTMin), _dTMax);
		
		if (newDeltaT <= _dT)
		{
			// if there is a drop in dt, just use the new value
			// (also do this if it is not changing b/c it is at the minimum)
			_dT = newDeltaT;
			
			// in this case, we also need to re-populate the saved steps
			// so that we don't use too-large values to raise the step again
			for (int i=0;i<10;i++)
				_oldStep[i] = _dT;
		}
		else
		{
			// otherwise gradually update the timestep in case of increase
			
			// make the new option the average of 10 previous steps
			for (int i=1;i<10;i++)
				_oldStep[i]=_oldStep[i-1];
			_oldStep[0] = newDeltaT;
			newDeltaT = utils.ExtraMath.mean(_oldStep);

			// again make sure the step isn't too small or too large
			_dT = _oldStep[0] = Math.min(Math.max(newDeltaT,_dTMin),_dTMax);
		}
		
		// make the step into a nicer number
		_dT = ((int)Math.floor(_dT/(_dTMin/10)))*(_dTMin/10);
		
		LogFile.writeLog("TimeStep "+_dT+" ("+tOpt+")");
	}

	/**
	 * \brief Resets the simulation timestep and counters to zero
	 * 
	 * Sets simulation time step and counters to zero.
	 */
	public static void reset() 
	{
		_now = 0;
		_nIter = 0;
	}

	/**
	 * \brief Return the time spent since the beginning of the simulation (in hours)
	 * 
	 * Utility method that returns the time spent since the beginning of the simulation (in hours)
	 * 
	 * @return Double value that specifies the time spent since the beginning of the simulation (in hours)
	 */
	public static double getCurrentTime() 
	{
		return _now;
	}

	/**
	 * \brief Returns the current iteration of the simulation
	 * 
	 * Returns the current iteration of the simulation
	 * 
	 * @return	 Integer value noting the current iteration of the simulation
	 */
	public static int getCurrentIter() 
	{
		return _nIter;
	}


	/**
	 * \brief Returns the currrent timestep, in hours
	 * 
	 * Returns the currrent timestep, in hours
	 * 
	 * @return Double value stating the current simulation timestep
	 */
	public static double getCurrentTimeStep()
	{
		return _dT;
	}

	/**
	 * \brief Sets the current timestep to a new double value. Used in adaptive timestep scenarios
	 * 
	 * Sets the current timestep to a new double value. Used in adaptive timestep scenarios
	 * 
	 * @param dt	The double value to assign to the simulation timestep
	 */
	public static void setCurrentTimeStep(double dt)
	{
		_dT = dt;
	}

	/**
	 * \brief Set the state of the simulation timer from information in an XML file
	 * 
	 * Set the state of the simulation timer from information in an XML file
	 * 
	 * @param infoFile	Name of the XML file containing the sim timer information
	 */
	public void setTimerState(String infoFile)
	{
		XMLParser fileRoot = new XMLParser(infoFile);
		XMLParser localRoot = new XMLParser(fileRoot.getChildElement("simulation"));
		
		_nIter = Integer.parseInt(localRoot.getAttribute("iterate"));
		_now = Double.parseDouble(localRoot.getAttribute("time"));
	}

	/**
	 * \brief Determine whether or not the simulation has finished
	 * 
	 * Determine whether or not the simulation has finished
	 * 
	 * @return	Boolean stating whether or not the simulation has finished
	 */
	public static boolean simIsFinished()
	{
		return (_now >= _endOfSimulation);
	}
	
	/**
	 * \brief Returns the length of the simulation in hours
	 * 
	 * Returns the length of the simulation in hours
	 * 
	 * @return Length of the simulation in hours
	 */
	public double returnSimulationLength()
	{
		return _endOfSimulation;
	}

	/**
	 * \brief Determines if a set time (such as an agent birthday) is within the next timestep
	 * 
	 * Determines if a set time (such as an agent birthday) is within the next timestep. Used to determine when agents enter the simulation
	 * 
	 * @param aDate	The timestep at which an event is set to happen
	 * @return	Boolean stating whether this event should occur in the next timestep
	 */
	public static boolean isDuringNextStep(double aDate) 
	{
		return (aDate >= _now) && (aDate < _now+_dT);
	}
	
	/**
	 * \brief Determines if the simulation time is within a set range of hours
	 * 
	 * Determines if the simulation time is within a set range of hours. Useful for determining whether an event can occur in a set 
	 * window (i.e. cell entry)
	 * 
	 * @param start	Hour at which event should start to occur
	 * @param end	Hour at which event should cease to occur
	 * @return	Boolean stating whether this event should occur in the next timestep
	 */
	public static boolean isCurrentTimeStepInSetInRange(double start, double end)
	{
		return (start <= _now) && (end >= _now);
	}
}
