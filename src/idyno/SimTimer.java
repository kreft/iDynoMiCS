/**
 * \package idyno
 * \brief Package of classes used to launch iDynomics.
 * 
 * Package of classes used to launch and iDynoMiCS simulation, and to update
 * the package to the latest stable release. This package is part of iDynoMiCS
 * v1.2, governed by the CeCILL license under French law and abides by the
 * rules of distribution of free software. You can use, modify and/ or
 * redistribute iDynoMiCS under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package idyno;

import simulator.World;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Class to create and keep track of the timestep and simulation time
 * course.
 *
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 */
public class SimTimer
{
	/**
	 * Number of simulation iterations performed.
	 */
	private static int _nIter;
	
	/**
	 * Time step of the simulation, in hours. Read in from the timeStepIni
	 * parameter of the XML protocol file.
	 */
	private static Double _dT;
	
	/**
	 * As _dT changes throughout the simulation, statedTimeStep holds the value
	 * that was originally read into _dT.
	 */
	public Double _statedTimeStep;
	
	/**
	 * Boolean flag that notes whether adaptive time steps are being employed.
	 * Read in from adaptive parameter in XML protocol file.
	 */
	private static Boolean isAdaptive;
	
	/**
	 * Maximum time step value (in hours) that can be used by the simulator,
	 * when adaptive time steps are employed. Read in from timeStepMax
	 * parameter in XML protocol file.
	 */
	private static Double _dTMax;
	
	/**
	 * Minimum time step value (in hours) that can be used by the simulator,
	 * when adaptive time steps are employed.  Read in from timeStepMin
	 * parameter in XML protocol file.
	 */
	private static Double _dTMin;
	
	
	private static Double[] _oldStep;
	
	/**
	 * Simulation time - in hours.
	 */
	private static Double _now;
	
	/**
	 * Amount of time for which the simulation should run (in hours).
	 */
	private static Double _endOfSimulation;
	
	
	/*************************************************************************
	 * CLASS METHODS 
	 *************************************************************************/
	
	
	
	/**
	 * \brief Parses simulation time step information from the XML protocol file
	 * 
	 * The timeStep markup within the protocol file specified all parameters
 	 * for the time-stepping of the simulation. This action is controlled by an
 	 * object of this SimTimer class. There is a choice of using an adaptive or
 	 * set timestep, specified by setting the adaptive parameter to true or
 	 * false respectively. If an adaptive timestep is used, three timeStep
 	 * parameters are set to control the initial, minimum, and maximum values
 	 * the timestep may take. The simulation then determines during runtime the
 	 * timestep that should be used. On the other hand, if this is set to
 	 * false, the simulation runs at one timestep - that set in the parameter 
 	 * timeStepIni. A further parameter is read in, endOfSimulation, that
  	 * specifies when the simulation should end.
	 * 
	 * @param localRoot	The Simulator mark-up within the XML file containing
 	 * the timestep parameters.
	 */
	public SimTimer(XMLParser localRoot) 
	{
		XMLParser parser = new XMLParser(localRoot.getChildElement("timeStep"));

		// Set all counters to zero.
		reset();

		// Initialize the timer.
		_endOfSimulation = parser.getParamTime("endOfSimulation");
		_statedTimeStep = _dT = parser.getParamTime("timeStepIni");

		// Now determine if adaptive timesteps are being used.
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
	 * \brief Computes the current simulation time step and increments
	 * the simulation clock.
	 * 
	 * First method called by the simulation step method.
	 */
	public static void applyTimeStep() 
	{
		_now += _dT;
		_nIter++;
	}

	/**
	 * \brief If adaptive timesteps are enabled, updates the timestep for a
	 * given simulation world.
	 * 
	 * @param aWorld The simulation world being modelled under an adaptive timestep.
	 */
	public static void updateTimeStep(World aWorld)
	{
		if ( !isAdaptive )
			return;
		
		Double tOpt, newDeltaT;
		
		tOpt = aWorld.getBulkTimeConstraint();
		
		if ( ! Double.isFinite(tOpt) )
			return;
		
		// Constrain value between min and max limits.
		newDeltaT = Math.min(Math.max(tOpt, _dTMin), _dTMax);
		
		// If dT needs to rise, increase it gradually. Otherwise, it must be
		// staying the same (no change) or instantly dropping to a new value.
		if (newDeltaT > _dT)
		{
			// Make the new option the average of 10 previous steps.
  			for (int i=1;i<10;i++)
  				_oldStep[i]=_oldStep[i-1];
  			_oldStep[0] = newDeltaT;
  			newDeltaT = ExtraMath.mean(_oldStep);
  		// Again make sure the step isn't too small or too large.
			_dT = _oldStep[0] = Math.min(Math.max(newDeltaT, _dTMin), _dTMax);
		}
		else
		{
			_dT = newDeltaT;
			
			// In this case, we also need to re-populate the saved steps so
			// that we don't use too-large values to raise the step again.
			for (int i = 0; i < 10; i++)
				_oldStep[i] = _dT;
		}
		
		// Make the step into a nicer number.
		_dT = Math.floor(10.0*_dT/_dTMin) * _dTMin * 0.1;
		
		LogFile.writeLog("TimeStep "+_dT+" ("+tOpt+")");
	}
	
	/**
	 * \brief Resets the simulation timestep and counters to zero.
	 */
	public static void reset() 
	{
		_now = 0.0;
		_nIter = 0;
	}

	/**
	 * \brief Return the time spent since the beginning of the simulation (in
	 * hours).
	 * 
	 * @return Double value that specifies the time spent since the beginning
	 * of the simulation (in hours).
	 */
	public static Double getCurrentTime() 
	{
		return _now;
	}

	/**
	 * \brief Returns the current iteration of the simulation.
	 * 
	 * @return	 Integer value noting the current iteration of the simulation.
	 */
	public static int getCurrentIter() 
	{
		return _nIter;
	}
	
	/**
	 * \brief Returns the current timestep, in hours.
	 * 
	 * @return Double value stating the current simulation timestep.
	 */
	public static Double getCurrentTimeStep()
	{
		return _dT;
	}

	/**
	 * \brief Sets the current timestep to a new double value.
	 * 
	 * Used in adaptive timestep scenarios.
	 * 
	 * @param dt	The double value to assign to the simulation timestep.
	 */
	public static void setCurrentTimeStep(Double dt)
	{
		_dT = dt;
	}
	
	/**
	 * \brief Set the state of the simulation timer from information in an XML file.
	 * 
	 * @param infoFile	Name of the XML file containing the sim timer information.
	 */
	public void setTimerState(String infoFile)
	{
		XMLParser fileRoot = new XMLParser(infoFile);
		XMLParser localRoot = new XMLParser(fileRoot.getChildElement("simulation"));
		_nIter = localRoot.getAttributeInt("iterate");
		_now = localRoot.getAttributeDbl("time");
	}

	/**
	 * \brief Determine whether or not the simulation has finished.
	 * 
	 * @return	Boolean stating whether or not the simulation has finished.
	 */
	public static Boolean simIsFinished()
	{
		return (_now >= _endOfSimulation);
	}
	
	/**
	 * \brief Returns the length of the simulation in hours.
	 * 
	 * @return Length of the simulation in hours.
	 */
	public Double returnSimulationLength()
	{
		return _endOfSimulation;
	}

	/**
	 * \brief Determines if a set time (such as an agent birthday) is within the next timestep.
	 * 
	 * Used to determine when agents enter the simulation.
	 * 
	 * @param aDate	The timestep at which an event is set to happen.
	 * @return	Boolean stating whether this event should occur in the next timestep.
	 */
	public static Boolean isDuringNextStep(Double aDate) 
	{
		return (aDate >= _now) && (aDate < _now+_dT);
	}
	
	/**
	 * \brief Determines if the simulation time is within a set range of hours.
	 * 
	 * Useful for determining whether an event can occur in a set 
	 * window (i.e. cell entry).
	 * 
	 * @param start	Hour at which event should start to occur.
	 * @param end	Hour at which event should cease to occur.
	 * @return	Boolean stating whether this event should occur in the next timestep.
	 */
	public static Boolean isCurrentTimeStepInSetInRange(Double start, Double end)
	{
		return (start <= _now) && (end >= _now);
	}
}
