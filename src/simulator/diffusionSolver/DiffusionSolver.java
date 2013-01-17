/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * DiffusionSolver is an abstract class used as parent for all diffusion_solvers 
 * you could define
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */

package simulator.diffusionSolver;

import java.util.*;

import org.jdom.Element;

import simulator.*;
import simulator.geometry.Domain;
import simulator.reaction.Reaction;

import utils.LogFile;
import utils.XMLParser;

public abstract class DiffusionSolver {

	public String                 solverName;
	public int                    solverIndex;
	public Simulator              mySim;
	public Domain    myDomain;

	// the list of pathways (reactions) to solve
	protected ArrayList<Reaction> _reactions   = new ArrayList<Reaction>();

	// List of solute grids created in simulator (just a reference)
	protected SoluteGrid[]        _soluteList;
	// List of index of solutes REALLY used by this solver
	protected ArrayList<Integer>  _soluteIndex = new ArrayList<Integer>();

	protected double              internTimeStep, minimalTimeStep;
	protected double              internalIteration;
	protected double              maxIteration = 1;

	protected Boolean			_active = false;

	/* ____________________ CREATION AND INITIALISATION _________________ */

	/**
	 * Initialisation procedure
	 * @param aSim
	 * @param xmlRoot
	 */
	public void init(Simulator aSim, XMLParser xmlRoot) {
		String reacName;

		mySim = aSim;
		solverName = xmlRoot.getAttribute("name");

		myDomain = aSim.world.getDomain(xmlRoot.getAttribute("domain"));

		_active = xmlRoot.getParamBool("active");

		// Reference all the solutes declared in this system
		_soluteList = aSim.soluteList;

		// Now add the reactions and list the solutes they modify
		LinkedList<Element> ReactionList = xmlRoot.buildSetMarkUp("reaction");
		for (Element aReactionMarkUp : ReactionList) {
			reacName = aReactionMarkUp.getAttributeValue("name");
			addReactionWithSolutes(aSim.getReaction(reacName));
		}
	}

	/**
	 * Sets reference to a biochemical pathway this solver has to deal with.
	 * References to the solutes and agents of the diffusion/reaction-system are
	 * provided by the pathways.
	 * @param aReaction : the pathway to add to this solver
	 */
	protected void addReactionWithSolutes(Reaction aReaction) {
		int aSoluteIndex;

		// assign the reactions
		_reactions.add(aReaction);

		// Collect references to SoluteGrids from the pathway and store them
		for (String aSoluteName : aReaction.declareSolutes()) {
			aSoluteIndex = mySim.getSoluteIndex(aSoluteName);
			if (!_soluteIndex.contains(aSoluteIndex)) {
				_soluteIndex.add(aSoluteIndex);
			}
		}
	}

	public void register() {
		try{
			solverIndex = mySim.getSolverIndex(solverName);
			mySim.solverList[solverIndex] = this;
		}catch(Exception e){LogFile.writeLog("Error in DiffusionSolver.register()");}
	}

	/**
	 * Small routine to use if you have only one solver instead to add one by
	 * one all pathways
	 */
	public void addAllReactions() {
		for (int i = 0; i<mySim.reactionList.length; i++) {
			addReactionWithSolutes(mySim.reactionList[i]);
		}
	}

	public Boolean isActive() {
		return _active;
	}

	/* ____________________ CREATION AND INITIALISATION _________________ */
	public void initAndSolve(){
		if (isActive()) {
			initializeConcentrationFields();
			solveDiffusionReaction();
		}
	}

	/**
	 * Initialize the diffusion-reaction-system according to the solver. Creates
	 * and initializes internal data structure for solving. Called at each
	 * simulation step
	 */
	public abstract void initializeConcentrationFields();

	/**
	 * Performs the solving algorithm on the diffusion reaction system. If
	 * needed, the time step is provided by the SimTimer.
	 */
	public abstract void solveDiffusionReaction();

}
