/**
 * \package diffusionSolver
 * \brief Package of classes used to capture the diffusion solvers that can be defined in the protocol file
 * 
 * Package of classes used to capture the diffusion solvers that can be defined in the protocol file. Solvers are used to compute 
 * the steady-state solute profile within the computational domains. This package is part of iDynoMiCS v1.2, governed by the CeCILL 
 * license under French law and abides by the rules of distribution of free software. You can use, modify and/ or redistribute iDynoMiCS 
 * under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.diffusionSolver;

import java.util.*;

import org.jdom.Element;

import simulator.*;
import simulator.geometry.Domain;
import simulator.reaction.Reaction;

import utils.LogFile;
import utils.XMLParser;

/**
 * 
 * \brief An abstract class used as a parent for all diffusion solvers that could be defined 
 * 
 * This class is used as a parent for all diffusion solvers that can be specified in the XML protocol file. This processes the SOLVER 
 * mark-up within the protocol file. Solvers are used to compute the steady-state solute profile within the computational domains.
 * 
 * This class is a component class of iDynoMiCS, released under the CECIL license. Please see www.idynomics.bham.ac.uk for more information
 *
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public abstract class DiffusionSolver 
{

	/**
	 * A name assigned to this solver. Specified in the XML protocol file
	 */
	public String                 solverName;
	
	/**
	 * The position of this solver in the simulation dictionary
	 */
	public int                    solverIndex;
	
	/**
	 * The simulation object that this solver is associated with
	 */
	public Simulator              mySim;
	
	/**
	 * The computational domain that this solver is associated with. Specified in the XML protocol file
	 */
	public Domain    myDomain;
	
	/**
	 * List of pathways (reactions) that this solver is associated with
	 */
	protected ArrayList<Reaction> _reactions   = new ArrayList<Reaction>();

	/**
	 * Local copy of the array of solute grids - one for each solute specified in the simulation protocol file. Taken from simulator object
	 */
	protected SoluteGrid[]        _soluteList;
	
	/**
	 * List of solutes that are used by THIS solver
	 */
	protected ArrayList<Integer>  _soluteIndex = new ArrayList<Integer>();

	protected double              internTimeStep;

	protected double              internalIteration;
	
	protected double              maxIteration = 1;

	/**
	 * Boolean flag that determines whether this solver will actually be used. Specified in XML protocol file
	 */
	protected Boolean			_active = false;
	
	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/

	/**
	 * \brief Initialisation procedure for each diffusion solver specified in the XML protocol file
	 * 
	 * This method takes a solver specification from the XML file (as a set of XML tags) and initialises a solver object
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param xmlRoot	The XML object containing the definition of one solver in the protocol file
	 */
	public void init(Simulator aSim, XMLParser xmlRoot) 
	{
		String reacName;

		mySim = aSim;
		
		// Get the name of this solver from the XML tag
		solverName = xmlRoot.getAttribute("name");

		// Get the computational domain that this solver is associated with
		myDomain = aSim.world.getDomain(xmlRoot.getAttribute("domain"));

		// Determine whether this solver will actually be used
		_active = xmlRoot.getParamBool("active");

		// Reference all the solutes declared in this system
		_soluteList = aSim.soluteList;

		// Now for each solver, reactions are specified. Add these reactions and list the solutes that these modify
		// Get the reaction tags for this solver
		LinkedList<Element> ReactionList = xmlRoot.buildSetMarkUp("reaction");
		
		// Now iterate through these
		for (Element aReactionMarkUp : ReactionList) 
		{
			// Get the name of the reaction that is modified
			reacName = aReactionMarkUp.getAttributeValue("name");
			
			// Now set the reference to a biochemical pathway this solver has to deal with
			addReactionWithSolutes(aSim.getReaction(reacName));
		}
	}

	/**
	 * \brief Sets reference to a biochemical pathway this solver has to deal with.
	 * 
	 * Sets reference to a biochemical pathway this solver has to deal with. References to the solutes and agents of the 
	 * diffusion/reaction-system are provided by the pathways.
	 * 
	 * @param aReaction : the pathway to add to this solver
	 */
	protected void addReactionWithSolutes(Reaction aReaction) 
	{
		// Used to store a reference to a solute in the simulation dictionary
		int aSoluteIndex;

		// assign the reactions
		_reactions.add(aReaction);

		// Collect references to SoluteGrids from the pathway and store them
		for (String aSoluteName : aReaction.declareSolutes()) 
		{
			// Get the index for this solute
			aSoluteIndex = mySim.getSoluteIndex(aSoluteName);
			
			// Add this to the list of solutes affected by this solver if not already present
			if (!_soluteIndex.contains(aSoluteIndex)) 
			{
				_soluteIndex.add(aSoluteIndex);
			}
		}
	}

	/**
	 * \brief Registers this solver in the simulation solver array for referencing later on
	 * 
	 * This method adds this solver object to the simulation variable solverList, an array of DiffusionSolver objects used in the 
	 * simulation
	 */
	public void register() 
	{
		try
		{
			// Get the position of this solver in the simulation dictionary
			solverIndex = mySim.getSolverIndex(solverName);
			
			// Add this to the simulation list of solvers
			mySim.solverList[solverIndex] = this;
		}
		catch(Exception e)
		{
			LogFile.writeLog("Error in DiffusionSolver.register()");
		}
	}

	/**
	 * \brief Small routine to use if you have only one solver instead to add one by one all pathways
	 * 
	 * Small routine to use if you have only one solver instead to add one by one all pathways
	 */
	public void addAllReactions() {
		for (int i = 0; i<mySim.reactionList.length; i++) {
			addReactionWithSolutes(mySim.reactionList[i]);
		}
	}

	/**
	 * \brief Determine if this solver is actually being used. Set in the protocol file
	 * 
	 * Determine if this solver is actually being used. Set in the protocol file
	 * 
	 * @return	Boolean stating whether this solver is active
	 */
	public Boolean isActive() 
	{
		return _active;
	}

	/**
	 * \brief Create the solver, initialise the concentration fields, and solve the diffusion reaction equations
	 * 
	 * Create the solver, initialise the concentration fields, and solve the diffusion reaction equations
	 */
	public void initAndSolve(){
		if (isActive()) {
			initializeConcentrationFields();
			solveDiffusionReaction();
		}
	}

	/**
	 * \brief	Initialize the diffusion-reaction-system according to the solver. 
	 * 
	 * Initialize the diffusion-reaction-system according to the solver. Creates and initializes internal data structure for solving. 
	 * Called at each simulation step
	 */
	public abstract void initializeConcentrationFields();

	/**
	 * \brief Performs the solving algorithm on the diffusion reaction system. If needed, the time step is provided by the SimTimer.
	 * 
	 * Performs the solving algorithm on the diffusion reaction system. If needed, the time step is provided by the SimTimer.
	 */
	public abstract void solveDiffusionReaction();

}
