/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation and their participation in relevant reactions
 * 
 * Package of utilities that create and manage agents in the simulation and their participation in relevant reactions. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.agent;
import java.util.List;
import java.util.Random;
import java.io.Serializable;
import java.awt.Color;

import org.jdom.Element;

import simulator.Simulator;
import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;
import utils.LogFile;
import utils.XMLParser;
import utils.ExtraMath;

/**
 * \brief Creates and manages the species that are to be included in an iDynoMiCS simulation
 * 
 * The Species class creates and manages the species that are to be included in an iDynoMiCS simulation. This includes the creation of 
 * the required number of agents, registering that this agent has been born, and reducing the count of active agents when the agent 
 * dies. The species each exist on their own grid, again initialised here. From version 1.2 of iDynoMiCS, this class also manages agent 
 * input for both one-time and self attachment scenarios.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 *
 */
public class Species implements Serializable 
{
	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long  serialVersionUID = 1L;

	/**
	 * Local copy of the simulation object that is creating the state specified in the protocol file
	 */
	public Simulator           currentSimulator;
	
	/**
	 * Name of the species associated with this object. Specified in the protocol file
	 */
	public String              speciesName;
	
	/**
	 * The number of species in the simulation
	 */
	public int                 speciesIndex;
	
	/**
	 * Colour used to distinguish this species in POV-Ray output images
	 */
	public Color               color;

	/**
	 *  Computational domain that this species is associated with
	 */
	public Domain              domain;
	
	/**
	 * Stores the position of agents that are swimming below the boundary until they meet the surface or the biofilm. Used in self-attach
	 * simulations
	 */
	public ContinuousVector 	swimmingAgentPosition;
	
	/**
	 * For self-attachment scenarios, the hour at which the injection of cells into the domain will start
	 */
	public double				cellInjectionStartHour;
	
	/**
	 * For self-attachment scenarios, the frequency (agents per hour) that attach to the substratum surface when entering from the  
	 * boundary layer, when cells are being injected into the domain
	 */
	public double 				injectionOnAttachmentFrequency;

	/**
	 * For self-attachment scenarios, the hour at which the injection of cells into the domain will stop
	 */
	public double 				cellInjectionEndHour;
	
	/**
	 * For self-attachment scenarios, the frequency (agents per hour) that attach to the substratum surface when entering from the  
	 * boundary layer, when cells are being injected into the domain
	 */
	public double 				injectionOffAttachmentFrequency;
	
	/**
	 * Used in self-attachment cases where the number of agents specified in the protocol file leaves a remainder per timestep. When this 
	 * remainder reaches 1, a new cell is introduced to the simulation
	 */
	public double				newAgentCounter;
	
	/**
	 * For self-attaching species, this holds the XY angle the cell is moving. This is global as this may change due to bouncing of 
	 * boundaries
	 */
	public double				angleOfMovingAgentXY;
	
	/**
	 * For self-attaching species, this holds the XZ angle the cell is moving (if 3D). This is global as this may change due to bouncing of 
	 * boundaries
	 */
	public double				angleOfMovingAgentXZ;

	/**
	 * Specialised agent from which objects of this species are created
	 */
	protected SpecialisedAgent _progenitor;
	
	/**
	 * Count of the population of this particular species in the simulation
	 */
	protected int              _population      = 0;

	/**
	 * \brief Creates a species object from that specified in the XML protocol file
	 * 
	 * This method takes an object specified in the SPECIES mark-up and creates a simulation object for that species
	 * 
	 * @param aSimulator	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aSpRoot	A Species mark-up within the specified protocol file
	 */
	public Species(Simulator aSimulator, XMLParser aSpRoot, XMLParser speciesDefaults) 
	{
		// Name of the species as specified in the protocol file
		speciesName = aSpRoot.getAttribute("name");
		
		// colour is used to distinguish agents in POV-Ray output images - read this from the protocol file
		String colorName = aSpRoot.getParam("color");
		
		// Set the colour to white if not specified
		if (colorName==null) colorName = "white";
		// Translate this text string into a colour
		color = utils.UnitConverter.getColor(colorName);

		// Register this species
		// KA 8/3/13 - is this correct - this is storing the number of species in the simulation
		speciesIndex = aSimulator.speciesList.size();
		
		// Take a local copy of the current simulation
		currentSimulator = aSimulator;
		
		// KA May 2013
		// If this is self-attach, useful to store the injection period parameters here, so can reference these later
		cellInjectionStartHour = aSpRoot.getParamDbl("cellInjectionStartHour");
		injectionOnAttachmentFrequency = aSpRoot.getParamDbl("injectionOnAttachmentFrequency");
		cellInjectionEndHour = aSpRoot.getParamDbl("cellInjectionEndHour");
		injectionOffAttachmentFrequency = aSpRoot.getParamDbl("injectionOffAttachmentFrequency");

		// Create the progenitor and tune its speciesParam object
		_progenitor = (SpecialisedAgent) aSpRoot.instanceCreator("simulator.agent.zoo");
		// Get parameters for this progenitor object from the protocol file if present

		_progenitor.getSpeciesParam().init(aSimulator, aSpRoot, speciesDefaults);
		
		_progenitor.setSpecies(this);

		// Set the computational domain this species is associated with
		// KA Aug 13 - changed as this may be a default
		domain = aSimulator.world.getDomain(_progenitor.getSpeciesParam().getSpeciesParameterString("computationDomain", aSpRoot, speciesDefaults));
	}

	/**
	 * \brief Create a species from a set progenitor
	 * 
	 * Create a species from a set progenitor
	 * 
	 * @param aProgenitor	Progenitor from which the species is being created
	 */
	public Species(SpecialisedAgent aProgenitor) 
	{
		_progenitor = aProgenitor;
		aProgenitor.setSpecies(this);
	}

	/**
	 * \brief Registers the created species to the simulation speciesManager
	 * 
	 * Registers the created species to the simulation speciesManager
	 * 
	 * @param aSimulator	The simulator object being used to create the conditions specified in the protocol file
	 * @param aSpRoot	XML markup for the species being created. Taken from the protocol file
	 */
	public void register(Simulator aSimulator, XMLParser aSpRoot) {
		currentSimulator = aSimulator;
		speciesName = aSpRoot.getAttribute("name");

		speciesIndex = aSimulator.speciesList.size();
		aSimulator.speciesList.add(this);

		domain = aSimulator.world.getDomain(aSpRoot.getAttribute("computationDomain"));
	}

	/**
	 * \brief Used in 'one-time' attachment scenarios, where clones of the progenitor are created in the birth area of the substratum
	 * 
	 * Used in 'one-time' attachment scenarios, where clones of the progenitor are created in the birth area of the substratum
	 * 
	 * @param spRoot	The XML mark-up group for a particular species being created
	 */
	public void createPop(XMLParser spRoot) 
	{
		double howMany = spRoot.getAttributeDbl("number");
	
		// Define the birth area - this is taken from the coordinates tags in the protocol file
		// (Nov13) OR if an initial area is not declared, this is the whole Y and Z of the domain with a height of 1.
		ContinuousVector[] _initArea = defineSquareArea(spRoot);
		String logStatement;
		
		if(Double.isNaN(howMany))
		{
			// User has declared a set number to be added over a specified area (for 3D) or length (for 2D)
			// We assume the user has declared the number of cells they want in millimetres in both cases
						
			// Now we need the area of the initial area / whole domain (depending how this has been declared) to calculate the number
			// of agents required. Make sure this is calculated in mm not microns
			if(domain.is3D) 
			{
				double cellsPerArea = spRoot.getAttributeDbl("cellsperMM2");
				double birthRegionArea = ((_initArea[1].y-_initArea[0].y)/1000)
										*((_initArea[1].z-_initArea[0].z)/1000);
				howMany = Math.floor(cellsPerArea * birthRegionArea);
				
				logStatement = "by specifying "+cellsPerArea+" cells per mm2";
			}
			else
			{
				
				double cellsPerLength = spRoot.getAttributeDbl("cellsperMM");
				double birthRegionLength = ((_initArea[1].y-_initArea[0].y)/1000);
				howMany = Math.floor(cellsPerLength * birthRegionLength);
				logStatement = "by specifying "+cellsPerLength+" cells per mm";
			}			
		}
		else
		{
			
			logStatement = "by specifying "+howMany+" cells for a stated initial area";
		}
		
		// Now we have a number and a birth area, we can add the cells
		if (howMany==0) return;

		// Create all the required agents
		ContinuousVector cc;
			
		for (int i = 0; i<howMany; i++) 
		{
			if (_progenitor instanceof LocatedAgent) 
			{
				cc = new ContinuousVector();
							
				if(!Simulator.isChemostat)
				{
					// Set coordinates within the birth area - randomly
					shuffleCoordinates(cc, _initArea);
					
				}
				
				// Create the agent at these coordinates
				((LocatedAgent) _progenitor).createNewAgent(cc);
				
			} 
			else 
			{
				_progenitor.createNewAgent();
			}
		}

		LogFile.writeLog(howMany+" agents of species "+speciesName+" for one-time attachment successfully created, "+logStatement);
	}
	
	
	
			
	/**
	 * \brief For self-attachment scenarios, initialises agents on the boundary layer rather than substrarum, and models their swim to the surface or biofilm
	 * 
	 * For self-attachment scenarios, the agents are initialised at the top of the boundary layer rather than on the substratum. These agents then 
	 * perform a 'run and tumble' motion until they either attach to the substratum or forming biofilm. This method captures this behaviour for 
	 * cells that are created for a time step. Once this swimming action has been performed, the agent is created at its final position. Note that 
	 * input of agents onto the boundary layer is decided by a parameter set in the protocol file, cellAttachmentFrequency, measured in hours. The 
	 * number of cells is adjusted to suit the global time step that is being used. Also note that this injection of cells can be for a set period 
	 * (specified in the protocol file as parameter cellInjectionPeriod), or can be stopped and started (modelling a 'settling' period) using 
	 * parameters cellInjectionOffPeriod and cellInjectionStopHour. This is explained in detail in the tutorial for version 1.2 of iDynoMiCS.
	 * 
	 * @param spRoot	The Species markup from the protocol file for one particular species being initialised
	 * @param numberAttachedInjectedAgents	The number of agents of this type that need to be created in this global timestep
	 */
	public void createBoundaryLayerPop(XMLParser spRoot, double numberAttachedInjectedAgents)
	{
		// Create all the required agents
		
		// Note that this continues UNTIL THE DESIRED NUMBER OF CELLS HAVE ATTACHED SOMEWHERE
		// Just out of interest, I've decided to keep a count of how many cells are required for this to happen
		int totalNumberOfInjectedAgents = 0;
		int agentsReturnedToBulk = 0;
		double requiredNumAttachedAgents = numberAttachedInjectedAgents;
		
		while(numberAttachedInjectedAgents > 0) 
		{
			totalNumberOfInjectedAgents++;
			
			if (_progenitor instanceof LocatedAgent) 
			{
				this.swimmingAgentPosition = new ContinuousVector();
				
				// Now to choose coordinates for this particular agent
				// boolean used to check that the coordinate chosen is ok
				boolean test = true;
				while (test) 
				{
					// This cell needs to take a random location in the Z and Y directions. The X will come from the position of the boundary layer 
					// on those axis. Generate these random numbers
					this.swimmingAgentPosition.y = ExtraMath.getUniRandDbl()*domain.length_Y;
					
					if(domain.is3D)
						this.swimmingAgentPosition.z = ExtraMath.getUniRandDbl()*domain.length_Z;
					else
						this.swimmingAgentPosition.z = 0;
				
					// Now to work out the X coordinate. This is based on where the top of the boundary layer is when this agent is created.
					// The top of the boundary layer is calculated in Domain at each step. Now the resolution differs (this is in nIxnJxnK rather than 
					// microns - so this will need to be converted accordingly
					// Formula to calculate this - get the value from the top of the boundary layer. Reduce by 1 (such that the micron value will be 
					// the top of the layer. Multiply by resolution of this simulation
					this.swimmingAgentPosition.x = ((domain._topOfBoundaryLayer
							[((int) Math.ceil((this.swimmingAgentPosition.y/domain._resolution)))+1]
							[((int) Math.ceil((this.swimmingAgentPosition.z/domain._resolution)))+1])-1)*domain._resolution;
							
					// Check this is ok
					test = !(domain.testCrossedBoundary(this.swimmingAgentPosition)==null);
				}
				
				// Now we can do the run and tumble motion of these cells
				
				int cellRunResult = performRunAndTumble(spRoot);
				
				
				// If increment the relevant counters, as these may be useful
				switch(cellRunResult)
				{
					case 1:	// Successfully Attached
						numberAttachedInjectedAgents--;						
						// Create the agent at these coordinates
						((LocatedAgent) _progenitor).createNewAgent(this.swimmingAgentPosition);
						break;
					case 2:
						agentsReturnedToBulk++;
						break;
				}
				
				
			}
			else 
			{
				_progenitor.createNewAgent();
			}
		}
		// Write the stats to the log file incase of interest
		LogFile.writeLog(requiredNumAttachedAgents+" agents of species "+speciesName+" for self-attachment successfully created");
		LogFile.writeLog(totalNumberOfInjectedAgents+" agents of species "+speciesName+" attempted to attach");
		LogFile.writeLog(agentsReturnedToBulk+" agents of species "+speciesName+" returned to the bulk");
		
	}
	
	/**
	 * \brief For self-attachment scenarios, performs the run and tumble motion of a swimming cell, and the required position checks and corrections
	 * 
	 * For self-attachment scenarios, the agents enter the system on top of the boundary layer and perform a run and tumble motion, swimming through 
	 * the biofilm. This method captures that motion. A random angle is chosen and the cell moves at a certain distance, set by cell tumble speed in 
	 * the protocol file. The agent may move back into the bulk, where we assume it will go back into the mix and thus we take no further action with 
	 * it. The aim is for that cell to swim around until it either hits the substratum surface or biofilm surface, where it will attach. If the agent 
	 * does not attach in this distance, a new angle is chosen (a 'tumble') and this is repeated. This method takes care of all checks to the coordinates 
	 * (such as a crossed boundary and hitting the surface/biofilm). An integer value is returned that states the fate of this move - the cell has 
	 * either not attached anywhere (0), and thus needs to tumble, the cell has attached (1), or the cell has returned to the bulk (2).
	 * 
	 * @author Kieran Alden
	 * @param spRoot	The Species markup from the protocol file for one particular species being initialised
	 * @return	 Integer value noting the fate of this move (0 no attachment, 1 attachment, 2 returned to bulk)
	 */
	public int performRunAndTumble(XMLParser spRoot)
	{
		// To clarify how Jan has designed this to work:
		// The moving cell has a speed, measured in microns a second, and this is specified in the protocol file.
		// You would thus think that the cell would move at a speed that this related to this and the global timestep.
		// However, it has been decided that this is not to be the case, and the cell will run and tumble continuously this timestep 
		// until it has attached to a surface or to the biofilm.
		// Thus this speed is only used to calculate how far the cell will move before tumbling, and the speed at which the cell moves
		// has no relation to the global timestep (and thus no relation to real time)
		
		// Now the distance the cell will move in a particular direction is calculated from the tumble interval (specified in the protocol file)
		// and the cell run speed (also in the protocol file). The cell will move away at a random angle for that distance, and will then 
		// tumble (changing the angle). This continues until this cell has either returned into the bulk or has attached to the surface
		// or biofilm
		// Eventually we will look to change this tumble interval, as this can be used to alter cell movement when under chemotaxis
		
		double distanceEachTumble = ((LocatedAgent) _progenitor).getSpeciesParam().cellRunSpeed *
									((LocatedAgent) _progenitor).getSpeciesParam().tumbleInterval;
		
		// We will use a integer flag system to determine what has happened to the cell after a run
		// 0 - no attachment, must do another run
		// 1 - cell has attached and thus movement is over
		// 2 - cell has returned into the bulk (and thus is discarded)
		int cellRunResult = 0;
		
		while(cellRunResult==0)
		{
			// Now work out the new position of the cell based on the distance covered this tumble
			// Firstly calculate that direction - generate a random angle between 0 and 360
			this.setAgentAngleOfMovement();
				
			// Now we're going to break the move down into chunks such that collision or move over the boundary can be detected
			// The size of the chunk is going to be the resolution of the agent grid
				
			double distanceMovedThisTumble = distanceEachTumble;
			while(distanceMovedThisTumble > 0 && cellRunResult==0)
			{
				// If the cell has less distance to move remaining that the resolution, should only move that much
				double distanceAgentMoves;
				
				// KA 03613 - move is radius of the cell + it's stickiness
				if(distanceMovedThisTumble > (((LocatedAgent) _progenitor).getSpeciesParam().divRadius 
				+ ((LocatedAgent) _progenitor).getSpeciesParam().stickinessAddition))
				{
					distanceAgentMoves = (((LocatedAgent) _progenitor).getSpeciesParam().divRadius 
							+ ((LocatedAgent) _progenitor).getSpeciesParam().stickinessAddition);
				}
				else
				{
					distanceAgentMoves = distanceMovedThisTumble;
				}
				
					
				calculateNewAgentPosition(distanceAgentMoves);
				
				// Now need to check for collision on any point of that path
				cellRunResult = checkAgentMove(distanceAgentMoves);
				
				// Subtract the distance moved in this mini move
				distanceMovedThisTumble = distanceMovedThisTumble - distanceAgentMoves;
			}
			
		}
		return cellRunResult;
	}
	
	/**
	 * \brief Set the angles an agent is to move in this tumble. If the domain is 3D, a move will involve calculation of
	 * 
	 *  Set the angles an agent is to move in this tumble. If the domain is 3D, a move will involve calculation of two angles. This is 
	 *  done using a random number generator 
	 */
	public void setAgentAngleOfMovement()
	{
		Random r = new Random();
		
		this.angleOfMovingAgentXY = 0 + 360 * r.nextDouble();
		
		if(domain.is3D)
			this.angleOfMovingAgentXZ = 0 + 360 * r.nextDouble();
		else
			this.angleOfMovingAgentXZ = 0;
	}
	
	/**
	 * \brief For self-attachment scenarios, calculates the new position of the agent based on XY and XZ angles and a distance to move
	 * 
	 * For self-attachment scenarios, calculates the new position of the agent based on XY and XZ angles and a distance to move. No return value 
	 * as the global swimmingAgentPosition is altered. The angles XY and XZ are also global parameters as these can also be altered by other 
	 * methods in the swimming agent position checks
	 * 
	 * @param distanceAgentMoves	Distance that the agent is to move
	 */
	public void calculateNewAgentPosition(double distanceAgentMoves)
	{
		// Now calculate where this angle would place the agent
		// Angle must be converted from degrees to radians for this to work
		this.swimmingAgentPosition.x = this.swimmingAgentPosition.x + distanceAgentMoves * 
				Math.cos(Math.toRadians(this.angleOfMovingAgentXY)) * Math.cos(Math.toRadians(this.angleOfMovingAgentXZ));
		this.swimmingAgentPosition.y = this.swimmingAgentPosition.y + distanceAgentMoves * Math.sin(Math.toRadians(this.angleOfMovingAgentXY));
		this.swimmingAgentPosition.z = this.swimmingAgentPosition.z + distanceAgentMoves * Math.sin(Math.toRadians(this.angleOfMovingAgentXZ));

	}
	
	/**
	 * \brief For self-attachment scenarios, determines whether a swimming agents move has taken it across a boundary, correcting the coordinate accordingly
	 * 
	 * For self-attachment scenarios, the agents are swimming through the domain, and we need to ensure that they perform the correct behaviour 
	 * when the boundary is met. This method checks whether an agents move has taken it over the boundary and applies the relevant correction (either 
	 * a bounce off the boundary or a reappearance on the other side). If the cell has hit the surface, the cell is deemed to have adhered to that 
	 * surface and a relevant x coordinate generated to note that this is the case. The top of the domain is dealt with differently, as this is 
	 * checked by the call to isNewCoordAboveBoundaryLayer, which determines if the agent has moved out of the boundary layer. If this is the case 
	 * we assume this cell to have returned to the bulk and do no further action with that cell. An integer is returned noting the fate of this 
	 * move - a 0 if the move is ok (after adjustment if required), a 1 if the agent has met the substratum and attached, and a 2 if the cell has 
	 * returned to the bulk 
	 * 
	 * @param distanceMoving	Distance the agent is moving (in microns) in this move
	 * @return	Integer noting the fate of this move (0 move ok (after adjustment if required), 1 if attached to surface, 2 if returned to bulk
	 */
	public int agentMoveBorderCheck(double distanceMoving)
	{
		// Simplest test to do first is to check if any boundaries have been crossed
		AllBC boundaryCrossed = domain.testCrossedBoundary(this.swimmingAgentPosition);
	
		// First is a simple test - has the cell returned to the bulk
		// If this is the case, we can just forget this and have another go
		// The cell will only have the capability to return to the bulk if the angle has it moving upwards or directly across
		// (i.e 0-90 and 270-360 degrees)
		if(!isNewCoordAboveBoundaryLayer())
		{
			// Now to see if the move takes the point outside any of the boundaries
		
			if(boundaryCrossed != null)
			{
				if(boundaryCrossed.getSideName().equals("y0z"))
				{
					// Detected that the move has crossed the substratum, thus the cell has hit the biofilm. 
					// A return of 1 indicates that this is the case
					// Hit the biofilm, so set the species coordinates as required
					
					// We may have hit the biofilm but the Y and Z coordinates in this generated move may still be negative (as they may have 
					// gone over another boundary. So before we set the final x, we should check Y and Z
					// So firstly, set the X position onto the surface
					this.swimmingAgentPosition.x = ExtraMath.getUniRandDbl();
					
					// Now set the final X position
					AllBC boundaryCrossedNewCheck = domain.testCrossedBoundarySelfAttach(this.swimmingAgentPosition);
					
					if(boundaryCrossedNewCheck != null)
					{
						if(boundaryCrossedNewCheck.getSideName().equals("xNz") || boundaryCrossedNewCheck.getSideName().equals("x0z"))
							this.correctCrossedLeftRightBoundaries(boundaryCrossedNewCheck, distanceMoving);
						else
							this.correctCrossedFrontBackBoundaries(boundaryCrossedNewCheck, distanceMoving);
						
					}
					
					
					return 1;
				}
				else if(boundaryCrossed.getSideName().equals("xNz") || boundaryCrossed.getSideName().equals("x0z"))
				{
					this.correctCrossedLeftRightBoundaries(boundaryCrossed, distanceMoving);
					
					return 0;
				}
				// Deal with 3D boundary too
				else if (boundaryCrossed.getSideName().equals("x0y") || boundaryCrossed.getSideName().equals("xNy"))
				{
					this.correctCrossedFrontBackBoundaries(boundaryCrossed, distanceMoving);
					
					return 0;
					
				}
				else
				{
					// This needs to be here so the function returns something. However this deals with the top boundary (yNz) and this
					// has already been dealt with by the crossed bulk method. So it is highly doubtful we will ever end up here
					 return 0;
				}
			}
			else
			{
				// No borders crossed, not back in the bulk
				return 0;
			}
		}
		else
		{
			// The cell has returned into the bulk, and thus this try is over. Return 2 so the process starts with a new cell
			return 2;
		}	
	}						
	
	/**
	 * \brief Used for self-attachment scenarios, to check the move of a swimming agent (crossed boundaries, hit surface, or hit biofilm)
	 * 
	 * In self-attachment scenarios, the cells move about a lot in each step. This method is employed to check whether the determined move is 
	 * suitable. Any move may take the agent across a boundary, which needs to be detected, and across the biofilm surface, which is what we are 
	 * aiming to hit. We also need to detect any contact with agents that are already part of a forming biofilm. This method calls others that 
	 * checks whether the generated move is suitable. The return is an integer value: 0 noting that the move is ok but there is no attachment 
	 * (to biofilm or surface), 1 that notes the agent has met the surface or biofilm, and 2 noting that the agent has returned to the bulk
	 * 
	 * @param distanceMoving	Distance that the swimming agent is moving in this particular move
	 * @return	Integer value denoting the result of this move (0 no attachment, 1 attachment, 2 return to bulk)
	 */
	public int checkAgentMove(double distanceMoving)
	{
		// Firstly, check if this move crosses any of the boundaries
		int checkBoundaries = agentMoveBorderCheck(distanceMoving);
		
		if(checkBoundaries == 1)
		{
			// Biofilm contact, no need to do any more
			return 1;
		}
		else if(checkBoundaries == 2)
		{
			// Returned to bulk 
			return 2;
		}
		else
		{
			// The agent has not crossed boundaries, yet has not yet stuck
			// Thus we need to determine if it is near a biofilm
			int bioFilmAttachment = this.checkForBioFilmContact(distanceMoving);
			
			return bioFilmAttachment;
		}
	}
	
	/**
	 * \brief Utility for determining how far a moving agent has moved past a boundary. Used in cyclic boundaries
	 * 
	 * Utility for determining how far a moving agent has moved past a boundary. Used in cyclic boundaries where if an agent has moved 
	 * past a boundary yet should appear that distance on the other side
	 * 
	 * @param newCoordinate	The new coordinate that is over the boundary
	 * @param lengthOfSide	The length of that side of the domain
	 * @param boundaryCrossed	The boundary that has been crossed in this agents move
	 * @return	Double noting the distance that this agent should reappear on the opposite side
	 */
	public double determineDistanceOverBoundary(double newCoordinate,double lengthOfSide,AllBC boundaryCrossed)
	{
		if(boundaryCrossed.getSideName().equals("xNy") || boundaryCrossed.getSideName().equals("xNz"))
			return newCoordinate - lengthOfSide;
		else
		{
			if(newCoordinate>=0)
				return newCoordinate - lengthOfSide;
			else
				return lengthOfSide + newCoordinate;
		}
	}
	
	/**
	 * \brief For self-attach scenarios, corrects 2D coordinates that have crossed the left of right boundary (x0z and xNz)
	 * 
	 * For self-attach scenarios, corrects 2D coordinates that have crossed the left of right boundary (x0z and xNz). The coordinate is 
	 * either placed on the opposite side for cyclic boundaries, or bounces off the boundary at the required angle
	 * 
	 * @param boundaryCrossed	The boundary that has been detected to have been crossed
	 * @param distanceMoving	The distance that the cell is moving in this step of its repositioning
	 */
	public void correctCrossedLeftRightBoundaries(AllBC boundaryCrossed, double distanceMoving)
	{
		// Cell has passed through the boundary on the left or right hand side
		// If is cyclic, the point is deemed to reappear at the opposite side
		if(boundaryCrossed.isCyclic())
		{
			// Calculate new point on the opposite side. 
			// This will be inside the grid at the amount that the move would have taken the point outside of the grid
			
			// Get the difference in J coordinate
			this.swimmingAgentPosition.y = determineDistanceOverBoundary(this.swimmingAgentPosition.y,(currentSimulator.agentGrid.get_nJ()*domain._resolution),boundaryCrossed);
		}
		else
		{
			// Can simply take the correct value of the new Y coordinate and reflect it the correct side of the boundary
			if(this.swimmingAgentPosition.y<0)
				this.swimmingAgentPosition.y = Math.abs(this.swimmingAgentPosition.y);
			else
				// gone over the other boundary, max size of the domain
				this.swimmingAgentPosition.y = ((domain._nJ*domain._resolution) - (this.swimmingAgentPosition.y - (domain._nJ*domain._resolution))-1);
			
			// Now we need to change the angle to note the change of direction
			this.angleOfMovingAgentXY = 360 - this.angleOfMovingAgentXY;
			
			// Worth checking the Z here if we're 3D - as we may have gone over a corner, and these may need reflecting in too
			if(domain.is3D)
			{
				if(this.swimmingAgentPosition.z<0)
					this.swimmingAgentPosition.z = Math.abs(this.swimmingAgentPosition.z);
				else
					// gone over the other boundary, max size of the domain
					this.swimmingAgentPosition.z = (domain._nK*domain._resolution) - (this.swimmingAgentPosition.z - (domain._nK*domain._resolution));
				
				// Note the change in the Z angle to note we've reflected
				this.angleOfMovingAgentXZ = 360 - this.angleOfMovingAgentXZ;
			}
		}
	}
	
	/**
	 * \brief For self-attach scenarios, corrects 3D coordinates that have crossed the front or back boundary (x0y and xNy)
	 * 
	 * For self-attach scenarios, corrects 3D coordinates that have crossed the front or back boundary (x0y and xNy). The coordinate is 
	 * either placed on the opposite side for cyclic boundaries, or bounces off the boundary at the required angle
	 * 
	 * @param boundaryCrossed	The boundary that has been detected to have been crossed
	 * @param distanceMoving	The distance that the cell is moving in this step of its repositioning
	 */
	public void correctCrossedFrontBackBoundaries(AllBC boundaryCrossed, double distanceMoving)
	{
		// Cell has passed through the boundary on the front or back of the grid
		// If is cyclic, the point is deemed to reappear at the opposite side
		if(boundaryCrossed.isCyclic())
		{
			// Calculate new point on the opposite side. 
			// This will be inside the grid at the amount that the move would have taken the point outside of the grid
			
			// Get the difference in K coordinate 
			this.swimmingAgentPosition.z = determineDistanceOverBoundary(this.swimmingAgentPosition.z,(currentSimulator.agentGrid.get_nK()*domain._resolution),boundaryCrossed);
		}
		else
		{
			// Can simply take the correct value of the new Y coordinate and reflect it the correct side of the boundary
			if(this.swimmingAgentPosition.z<0)
				this.swimmingAgentPosition.z = Math.abs(this.swimmingAgentPosition.z);
			else
				// gone over the other boundary, max size of the domain
			this.swimmingAgentPosition.z = (domain._nK*domain._resolution) - (this.swimmingAgentPosition.z - (domain._nK*domain._resolution));
						
			// Now we need to change the angle to note the change of direction
			this.angleOfMovingAgentXZ = 360 - this.angleOfMovingAgentXZ;
						
			// Worth checking the Y here - as we may have gone over a corner, and these may need reflecting in too
			
			if(this.swimmingAgentPosition.y<0)
				this.swimmingAgentPosition.y = Math.abs(this.swimmingAgentPosition.y);
			else
				// gone over the other boundary, max size of the domain
				this.swimmingAgentPosition.y = (domain._nJ*domain	._resolution) - (this.swimmingAgentPosition.y - (domain._nJ*domain._resolution));
						
			// Note the change in the Z angle to note we've reflected
			this.angleOfMovingAgentXY = 360 - this.angleOfMovingAgentXY;
			
		}
		
	}
	
	/**
	 * \brief Determine if a swimming agent has left the boundary layer and returned to the bulk
	 * 
	 * For self attachment, an agent starts at the top of the boundary layer and swims in a random direction until it attaches somewhere. 
	 * If however it returns to the bulk, we assume this cell does not return. This method checks the move to determine if the cell is 
	 * moving in a trajectory that has returned it to the bulk 
	 * 
	 * @return	Boolean stating whether the cell has returned to the bulk
	 */
	public boolean isNewCoordAboveBoundaryLayer()
	{
		boolean returnedToBulk = false;
		
		// Get the top of the bulk - note we may need to correct here. The method that calls this calculates new positions. However this is 
		// run before it is determined whether any boundaries have been crossed (as it would be silly to check all boundaries when there is a 
		// high chance, especially at the start, that the cell may return into the bulk. Thus, a high Y or Z value may be sent in here, which 
		// when referenced to the boundary layer array, will cause an error. So if the Y or Z values are higher than the size of the grid, 
		// we will use the size of the grid as the array reference
		if(this.swimmingAgentPosition.y > domain._nJ*domain._resolution)
			this.swimmingAgentPosition.y = (domain._nJ*domain._resolution)-1;
		
		if(this.swimmingAgentPosition.z > domain._nK*domain._resolution)
			this.swimmingAgentPosition.z = (domain._nK*domain._resolution)-1;
			
		// Now calculate the top of the boundary layer at the y and z coordinate
		double boundaryAtThisPoint = 
			domain._topOfBoundaryLayer[((((int)(this.swimmingAgentPosition.y/domain._resolution)))+1)]
								      [((((int)(this.swimmingAgentPosition.z/domain._resolution)))+1)]*domain._resolution;
		
		
		// Check if we are above it
		if(boundaryAtThisPoint < this.swimmingAgentPosition.x)
		{
			returnedToBulk = true;
		}
				
		return returnedToBulk;
	}
	
	/**
	 * \brief Determines whether a swimming agent has come into contact with the forming biofilm
	 * 
	 * For self-attachment scenarios, iDynoMiCS will have so far checked whether the swimming agent has crossed any of the boundaries, or hit the
	 * substratum. Now we need to determine if the agent has met any other agents that are on the biofilm surface. If this is the case, we assume 
	 * that the agent sticks and becomes part of the biofilm structure. If not, the run and tumble motion continues. An integer is returned noting 
	 * the result - 0 if there is no attachment to any part of the biofilm, 1 if the agent is either in contact with the biofilm or has hit the 
	 * substratum surface, and 2 if the agents move results in a return to the bulk
	 * 
	 * @author Kieran Alden
	 * @param distanceMoving	The distance (in microns) that the agent is swimming through the domain
	 * @return	Integer between 0 and 2, noting no attachment (0), contact with the biofilm or surface (1), or return to the bulk (2)
	 */
	public int checkForBioFilmContact(double distanceMoving)
	{
		int attachedToBioFilm = 0;
		
		// Get the index of the voxel that the agent resides in
		int index = currentSimulator.agentGrid.getIndexedPosition(this.swimmingAgentPosition);
		// Get the status of this voxel. This will be a 1 if considered part of the biofilm
		int voxelStatus = currentSimulator.agentGrid.getVoxelStatus(index);
		
		while(voxelStatus == 1 && attachedToBioFilm == 0)
		{
			// The move hits the biofilm
			// We did think this was sufficient - however the cell may have entered a grid space deemed to be part of the biofilm, 
			// yet still be a large distance from any other agent. Thus its move needs to continue
			
			// Lets break it into smaller chunks, and determine when the agent is in contact with another. We can do this using the
			// agent radius + stickiness as a distance
			
			double distanceSeekingAgent = ((LocatedAgent) _progenitor).getSpeciesParam().divRadius 
										+ ((LocatedAgent) _progenitor).getSpeciesParam().stickinessAddition;
			
			// Now go through each agent and determine if this is within that distance from an agent
			if(isAgentInContactWithAgentInBiofilm(index, distanceSeekingAgent))
			{
				// can return 1 as there is contact with an agent on the biofilm surface
				// Set the coordinates as this final position
				//cc.x = cc_For_Checking.x; cc.y = cc_For_Checking.y; cc.z = cc_For_Checking.z;
				attachedToBioFilm = 1;
			}
			else
			{
				// We need to shift the agent position, and then give the same test a try
				this.calculateNewAgentPosition(distanceSeekingAgent);
				
				// just check that this has not crossed any boundaries
				int boundaryCheck = this.agentMoveBorderCheck(distanceMoving);
				if(boundaryCheck == 1)
					// this move has taken the agent directly onto the surface, where we assume this attaches
					attachedToBioFilm = 1;
				else if(boundaryCheck == 2)
					// this move has taken the agent back into the bulk, so we can discard this agent
					attachedToBioFilm = 2;
				else
				{
					// the move is ok, but we're still not yet attached anywhere
					
					// Now this move may have taken the agent into a different voxel - so we need to get the voxel status
					// If this is the case, and the voxel status is again 1, then this routine continues. Else run and tumble should 
					// recommence (as dictated by the while loop)
					index = currentSimulator.agentGrid.getIndexedPosition(this.swimmingAgentPosition);
					voxelStatus = currentSimulator.agentGrid.getVoxelStatus(index);
				}
			}
		}
		
		// The while loop has either ended as we are in a voxel that is not in 'biofilm' status, we have attached, or we are back in 
		// the bulk. Whichever, set the coordinates to be the final check and return the relevant flag integer
		
		return attachedToBioFilm;
	}
	
	/**
	 * \brief Examines all the agents within a voxel of the agent-grid, to determine if a swimming agent is in contact with an agent in the biofilm
	 * 
	 * To determine if a swimming cell is near contact with the biofilm surface, the agent checks the status of the agent grid voxel it resides 
	 * within. If this is 1, this is noted as having agents within it that are within the biofilm structure. The agent now examines how close it 
	 * is to each of these agents. If within a certain distance (radius + stickiness constant), then the agent is deemed to have adhered to the 
	 * structure. Any overlap will be addressed during shoving analysis. If not, the agent will perform another move.
	 * 
	 * @author Kieran Alden
	 * @param gridIndex	The index of the agent grid to be checked
	 * @param distanceSeekingAgent	The distance within which two cells are deemed to be in contact
	 * @return	Boolean noting whether the agent is in contact with an agent in the biofilm
	 */
	public boolean isAgentInContactWithAgentInBiofilm(int gridIndex, double distanceSeekingAgent)
	{
		// Get the agents within that grid voxel
		LocatedGroup agentsInGrid = currentSimulator.agentGrid.returnGroupInVoxel(gridIndex);
		boolean contactMade = false;
		
		// Now iterate through each one
		for(int j=0;j<agentsInGrid.group.size() && !contactMade;j++)
		{
			LocatedAgent aLoc = agentsInGrid.group.get(j);
	
			// Now we're going to work out how far we are from any of the agents in the grid. If we're close enough, move done
			// shoving can then sort out distance between the two cells
			// If not, we'll do another move
	
			double distanceBetweenAgents = Math.sqrt( Math.pow((aLoc.getLocation().x - this.swimmingAgentPosition.x),2) +
					                                  Math.pow((aLoc.getLocation().y - this.swimmingAgentPosition.y),2) + 
				                                      Math.pow((aLoc.getLocation().z - this.swimmingAgentPosition.z),2));
	
			if(distanceBetweenAgents<=distanceSeekingAgent)
			{
				contactMade = true;
			}
		}
		
		return contactMade;
	}

	/**
	 * \brief Increases the population of this species when one agent is added
	 * 
	 * Increases the population of this species when one agent is added
	 */
	public void notifyBirth() {
		_population++;
	}

	/**
	 * \brief Reduces the population of this species when one agent is removed
	 * 
	 * Reduces the population of this species when one agent is removed
	 */
	public void notifyDeath() {
		_population--;
	}

	/**
	 * \brief Returns a clone of the progenitor of this species
	 * 
	 * Returns a clone of the progenitor of this species. Throws an exception if this clone cannot be found
	 * 
	 * @return a clone of the progenitor
	 * @throws CloneNotSupportedException
	 */
	public SpecialisedAgent sendNewAgent() throws CloneNotSupportedException {
		return _progenitor.sendNewAgent();
	}

	/**
	 * \brief Returns the population count of this species
	 * 
	 * Returns the population count of this species
	 * 
	 * @return	Integer value noting the population of this species
	 */
	public int getPopulation() {
		return _population;
	}

	/**
	 * \brief Return the progenitor of this species object
	 * 
	 * Return the progenitor of this species object
	 * 
	 * @return	The progenitor of this species object
	 */
	public SpecialisedAgent getProgenitor() 
	{
		return _progenitor;
	}

	/**
	 * \brief Return the parameters associated with this species (a SpeciesParam object)
	 * 
	 * Return the parameters associated with this species. These are contained within a SpeciesParam object
	 * 
	 * @return	SpeciesParam object associated with this Species
	 */
	public SpeciesParam getSpeciesParam() 
	{
		return _progenitor.getSpeciesParam();
	}

	/**
	 * \brief Return the species object for a given species name
	 * 
	 * Return the species object for a given species name
	 * 
	 * @param speciesName	Text string containing the name of this species
	 * @return	Species object named with the given text string
	 */
	public Species getSpecies(String speciesName) 
	{
		return currentSimulator.speciesList.get(currentSimulator.getSpeciesIndex(speciesName));
	}

	/**
	 * \brief Defines a region of the computation domain where a new species may be created, using restrictions in the protocol file
	 * 
	 * Defines a region of the computation domain where a new species may be created, using restrictions in the protocol file. These 
	 * restrictions for a particular species are specified in 'coordinates' tags. There should be two such tags where this is used - one
	 * to start the restriction and one to end it. Example of use: \<coordinates x="0" y="0" z="0"/\> \<coordinates x="1" y="264" z="0"/\>. 
	 * This method will read these in and create an array that represents this area
	 * 
	 * @param spRoot	The information within the 'initArea' tags of the protocol file
	 * @return	A continuous vector representing the area of the domain specified in these tags
	 */
	public ContinuousVector[] defineSquareArea(XMLParser spRoot) 
	{
		List<Element> area = spRoot.getChildren("coordinates");
		
		ContinuousVector[] initArea = new ContinuousVector[2];
		initArea[0] = new ContinuousVector();
		initArea[1] = new ContinuousVector();
		
		// KA NOV 13 - CHANGED THIS, AS WE'RE GOING TO LET THE USER NOT DECLARE AN INITIAL AREA IF THEY WANT THE CELLS SPREAD ACROSS
		// THE WHOLE DOMAIN. THUS THIS NEEDS CHECKING AND FIXING
		if(area.size()>0)
		{
			// First Coordinate Tag
			ContinuousVector cc1 = new ContinuousVector((Element) area.get(0));
			// Second Coordinate Tag
			ContinuousVector cc2 = new ContinuousVector((Element) area.get(1));
	
			// Set each point
			initArea[0].x = Math.min(cc1.x, cc2.x);
			initArea[0].y = Math.min(cc1.y, cc2.y);
			initArea[0].z = Math.min(cc1.z, cc2.z);
			initArea[1].x = Math.max(cc1.x, cc2.x);
			initArea[1].y = Math.max(cc1.y, cc2.y);
			initArea[1].z = Math.max(cc1.z, cc2.z);
	
		}
		else
		{
			// NO INITIAL AREA HAS BEEN DECLARED, USE THE WHOLE SUBSTRATUM. NOTE THAT THE X (HEIGHT) COORDINATE IS SET TO 1 SO THE
			// CELLS ARE PLACED NEAR THE SUBSTRATUM
			// Set each point
			initArea[0].x = 0;
			initArea[0].y = 0;
			initArea[0].z = 0;
			initArea[1].x = 1.0;
			initArea[1].y = domain.length_Y;
			initArea[1].z = domain.length_Z;
		}
		
		// In the case of 2D simulation, the agent's z-coordinate is 0.
		if (!domain.is3D) 
		{
			initArea[0].z = 0;
			initArea[1].z = 0;
		}
		
		return initArea;
	}

	/**
	 * \brief Select random coordinates for the new agent within a restricted birth area
	 * 
	 * Select random coordinates for the new agent within a restricted birth area. This restricted area is set within the protocol file, 
	 * and defined as a ContinuousVector by the method defineSquareArea
	 * 
	 * @param cc	ContinuousVector that will hold the coordinates of this agent
	 * @param area	Area within which these coordinates should be restricted
	 */
	public void shuffleCoordinates(ContinuousVector cc, ContinuousVector[] area) 
	{
		do 
		{
			cc.x = ExtraMath.getUniRandDbl(area[0].x, area[1].x);
			cc.y = ExtraMath.getUniRandDbl(area[0].y, area[1].y);
			cc.z = ExtraMath.getUniRandDbl(area[0].z, area[1].z);
		} while ( domain.testCrossedBoundary(cc) != null );
		 
	}
}
