
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 */



package simulator.reaction;
import java.io.Serializable;
import java.util.*;
import Jama.Matrix;
import simulator.*;
import simulator.agent.*;
import simulator.geometry.Bulk;
import simulator.geometry.DiscreteVector;
import simulator.geometry.boundaryConditions.AllBC;
import utils.LogFile;
import utils.ResultFile;
import utils.XMLParser;

@SuppressWarnings("serial")
public abstract class Reaction implements Serializable {

	/* ________________________ PROPERTIES ________________________________ */
	// List of process created
	public String                     reactionName;
	public int                        reactionIndex;
	public int                        _catalystIndex;

	// Agents hosting this process
	protected LinkedList<ActiveAgent> _guild = new LinkedList<ActiveAgent>();

	// List of solute grids needed by THIS process (however the size of the
	// array is defined by the total number of solutes)
	protected SoluteGrid[]            _soluteList;
	//sonia 21.09.09 I changed the visibility of _mySoluteIndex to public to access it from the Solver_chemostat class
	/**
	 * A list of solute indices: if that solute's yield in this reaction is non-zero then it is on this list.
	 * Hence, this is the list of solutes which are affected by this reaction. 
	 */
	public int[]                   _mySoluteIndex;
	/** Whether or not the catalytic particle is affected by the reaction */
	public boolean					autocatalytic;
	
	protected SoluteGrid              _reacGrid, _guildGrid;
	public double[]                   totalUptake;
	//public double                     globalReactionRate;

	// Temporary storage used during computation of reaction rates
	protected double                  _specRate;
	// Buffer vectors used to communicate with a solver. If a solute is
	// not concerned by the current reaction, the uptake will be zero
	protected double[]                _diffUptakeRate;
	protected double[]                _uptakeRate;



	// Parameters of the reaction : yield and kinetic params
	// If a solute is not concerned by the current reaction, the yield will be
	// zero
	protected double[]                _soluteYield;
	protected double[]                _kineticParam;
	protected double[]                _particleYield;

	//sonia: 21-05-09
	protected String[]             _particleNameYield ;

	//sonia 21-04-10
	//dilution rate from associated bulk
	public double Dil;

	// Temporary variables
	static int                        nSolute;

	/* _________________________ CONSTRUCTOR________________________________ */

	/**
	 * Initialization procedure ; based on XML file Kinetic parameters will be
	 * managed by the subdefined method in the children class
	 * @param aSim
	 * @param aReactionRoot
	 * @see Simulator.createReaction()
	 */
	public void init(Simulator aSim, XMLParser aReactionRoot) {
		reactionName = aReactionRoot.getAttribute("name");


		nSolute = aSim.soluteList.length;
		int nParticulate = aSim.particleDic.size();

		// Create a simple array of all solutes
		_soluteList = aSim.soluteList;
		_reacGrid = new SoluteGrid(_soluteList[0]);
		_reacGrid.gridName = reactionName+"-rate";

		_guildGrid = new SoluteGrid(_soluteList[0]);
		_guildGrid.gridName = reactionName+"-pop";

		// Initialize buffer arrays
		_uptakeRate = new double[nSolute];
		_diffUptakeRate = new double[nSolute];
		totalUptake = new double[nSolute];

		// Initialize arrays storing reaction parameters
		_soluteYield = new double[nSolute];
		_particleYield = new double[nParticulate];
		_particleNameYield = new String[nParticulate];

		for (AllBC aBC : _reacGrid.getDomain().getAllBoundaries()){
			if (aBC.hasBulk()){
				Bulk aBulk = aBC.getBulk();
				if(aBulk.getName().equals("chemostat")){
					Dil = aBulk._D;
				}
			}	
		}

		// Extract the yields for solutes and particulates from the XML file
		fillParameters(aSim, aReactionRoot);
	}

	/**
	 * Called by an ActiveAgent to populate their reaction parameters vectors
	 * Populate only the yield vector, the kinetic vector will be filled by the
	 * subdefined method in the child classes
	 * @param anActiveAgent : the agent calling the method
	 * @param aReactionRoot : the root describing the reaction hosted by the
	 * agent in the XML file
	 */
	public void initFromAgent(ActiveAgent anAgent, Simulator aSim, XMLParser xmlRoot) {
		double yield;

		// Populate yield for solutes __________________________________
		XMLParser parser = new XMLParser(xmlRoot.getChildElement("yield"));
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++) {
			yield = parser.getParamSuchDbl("solute", _soluteList[iSolute].getName());
			if (!Double.isNaN(yield)) {
				anAgent.soluteYield[reactionIndex][iSolute] = yield;
			}
		}

		// Populate yields for particles _________________________________
		for (int iParticle = 0; iParticle<aSim.particleDic.size(); iParticle++) {
			yield = parser.getParamSuchDbl("particle", aSim.particleDic.get(iParticle));
			if (!Double.isNaN(yield)) {
				anAgent.particleYield[reactionIndex][iParticle] = yield;
			}
		}
	}

	/**
	 * End of the initialization procedure
	 * @param aReactionName
	 * @param aSimulator
	 * @see Simulator.addReaction()
	 */
	// Called by Simulator.createReactions()
	public void register(String aReactionName, Simulator aSimulator) {
		reactionName = aReactionName;
		reactionIndex = aSimulator.getReactionIndex(aReactionName);
	}

	/**
	 * Called by a solver instance to know the solutes affected by this reaction
	 * @return a list of solutes affected by this reaction
	 * @see DiffusionSolver.addReaction()
	 */
	public LinkedList<String> declareSolutes() {
		LinkedList<String> affectedSolutes = new LinkedList<String>();
		for (int index : _mySoluteIndex) {

			affectedSolutes.add(_soluteList[index].getName());
		}

		return affectedSolutes;
	}

	/**
	 * Used during the initialization
	 */
	public void fillParameters(Simulator aSim, XMLParser xmlRoot) {

		// Who is mediating the reaction ?
		String catalystName = xmlRoot.getAttributeStr("catalyzedBy");
		_catalystIndex = aSim.particleDic.indexOf(catalystName);

		// Populate yield for solutes __________________________________
		double yield;
		int nSolute = 0;
		XMLParser parser = new XMLParser(xmlRoot.getChildElement("yield"));
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++) {
			//yield = parser.getParamSuchDbl("solute", _soluteList[iSolute].getName());
			yield = parser.getParamDbl( _soluteList[iSolute].getName());
			if (!Double.isNaN(yield)) {
				nSolute++;
				_soluteYield[iSolute] = yield;
			} else {
				_soluteYield[iSolute] = 0.0;
			}
			//LogFile.writeLog("Solute "+iSolute+" has yield "+_soluteYield[iSolute]);
		}
		
		// Populate yields for particle _________________________________
		String particleName;
		for (int iParticle = 0; iParticle<aSim.particleDic.size(); iParticle++) {
			//yield = parser.getParamSuchDbl("particle", aSim.particleDic.get(iParticle));
			yield = parser.getParamDbl(aSim.particleDic.get(iParticle));
			particleName = aSim.particleDic.get(iParticle);
			if (!Double.isNaN(yield)) {
				_particleYield[iParticle] = yield;
				_particleNameYield[iParticle] = particleName;
			}
		}
		
		// Determine whether the catalytic particle is affected by the reaction,
		// and hence whether the reaction is autocatalytic or not.
		// This is important in ActiveAgent.grow()
		if (_particleYield[_catalystIndex]==0){
			autocatalytic = false;
			LogFile.writeLog("Reaction "+reactionIndex+" is not autocatalyic");
		} else {
			autocatalytic = true;
			LogFile.writeLog("Reaction "+reactionIndex+" is autocatalyic");
		}

		int jSolute = 0;
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++)
			if (_soluteYield[iSolute]!=0) jSolute++;

		_mySoluteIndex = new int[jSolute];
		jSolute = 0;
		// For each solute, if that solute's yield in this reaction is non-zero,
		// add the solute's index to _mySoluteIndex
		//LogFile.writeLog("Reaction "+this.reactionIndex+" is affected by these solutes:");
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++) {
			//LogFile.writeLog("Yield of solute "+iSolute+" is "+_soluteYield[iSolute]);
			if (_soluteYield[iSolute]!=0) {
				_mySoluteIndex[jSolute] = iSolute;
				jSolute++;
				//LogFile.writeLog(" "+iSolute);
			}
		}


	}

	/* ______________ METHODS FOR REACTION MANAGEMENT_________________________ */

	/**
	 * Register an agent among the guild of this pathway
	 * @param anAgent
	 */
	public void addAgent(ActiveAgent anAgent) {
		_guild.add(anAgent);
	}

	/**
	 * Remove an agent among the guild of this pathway
	 * @param anAgent
	 */
	public void removeAgent(ActiveAgent anAgent) {
		_guild.remove(anAgent);
	}

	/**
	 * Compute the total mass of all members of the guild
	 * Doesn't seem to be called (11th Jan 2012)
	 * @return the reacting mass
	 */
	public double getReactingMass() {
		double totalMass = 0;
		for (ActiveAgent anAgent : _guild) {
			totalMass += anAgent.getParticleMass(_catalystIndex);
		}
		return totalMass;
	}

	public double getUptakeRate(int soluteIndex) {
		return _uptakeRate[soluteIndex];
	}

	/* ________________ COMMUNICATION WITH THE SOLVER _____________________ */

	/**
	 * Mass growth-rate (in gX.h-1)
	 */
	public abstract double computeMassGrowthRate(ActiveAgent anAgent);
	
	public abstract double computeSpecGrowthRate(ActiveAgent anAgent);

	/* Specific growth-rates are independent of agent mass ______________ */
	public abstract void computeSpecificGrowthRate(ActiveAgent anAgent);

	public abstract void computeSpecificGrowthRate(double[] s, ActiveAgent anAgent);

	public abstract void computeSpecificGrowthRate(double[] s);

	/**
	 * Add the contribution of this agent on the reaction grid and the diff
	 * reaction grid Catalyst quantity is expressed in CONCENTRATION
	 */
	public abstract void computeUptakeRate(double[] s, double mass, double h);

	public abstract void computeUptakeRate(double[] s, double mass, Matrix dFdY);

	/**
	 * Add the contribution of this agent on the reaction grid and the diff
	 * reaction grid
	 */
	public abstract void computeUptakeRate(double[] s, ActiveAgent anAgent);

	/**
	 * Compute reaction rate on each concerned solute grids Assumes same
	 * parameters for all the agents of a same guild
	 * @param concGrid : solute concentration
	 * @param reacGrid : contribution of the reaction to the solute
	 * concentration dynamics
	 * @param diffReacGrid : derivative of the previous grid
	 * @param biomassGrid : CONCENTRATION of the catalyst
	 * @see multigrid solver
	 */
	public void applyReaction(SpatialGrid[] concGrid, SpatialGrid[] reacGrid,
			SpatialGrid[] diffReacGrid, SpatialGrid biomassGrid) {

		nSolute = concGrid.length;
		double[] s = new double[nSolute];

		int _nI, _nJ, _nK;
		_nI = biomassGrid.getGridSizeI();
		_nJ = biomassGrid.getGridSizeJ();
		_nK = biomassGrid.getGridSizeK();
		//globalReactionRate = 0;


		for (int i = 1; i<_nI+1; i++) {
			for (int j = 1; j<_nJ+1; j++) {
				for (int k = 1; k<_nK+1; k++) {
					// If there is no biomass, go to the next grid element
					if (biomassGrid.grid[i][j][k]==0) continue;

					// Read local solute concentration
					for (int iGrid : _mySoluteIndex)
						s[iGrid] = concGrid[iGrid].grid[i][j][k];

					// First compute local uptake-rates in g.h-1
					//sonia:chemostat added a parameter to the computeUptakerate method
					computeUptakeRate(s, biomassGrid.grid[i][j][k], 0);
					//globalReactionRate += _specRate;

					// Now add them on the received grids
					for (int iGrid : _mySoluteIndex) {
						reacGrid[iGrid].grid[i][j][k] += _uptakeRate[iGrid];
						diffReacGrid[iGrid].grid[i][j][k] += _diffUptakeRate[iGrid];
						if (Double.isNaN(reacGrid[iGrid].grid[i][j][k])) 
							LogFile.writeLog("Warning: NaN generated in Reaction");
					}

				}
			}
		}
	}

	public abstract Matrix calcdMUdS(Matrix S, double biomass);
	public abstract Matrix calcdMUdT(Matrix S, double biomass);

	/**
	 * Does anything call this?
	 * @param concGrid
	 * @param reacGrid
	 */
	/*	public void applyReactionCA(SoluteGrid[] concGrid, SoluteGrid[] reacGrid,
	        SoluteGrid[] diffReacGrid) {
		double s[];
		double mass;

		if (true) {
			for (ActiveAgent anAgent : _guild) {
				s = readConcentrationSeen(anAgent, concGrid);
				mass = anAgent.getParticleMass(_catalystIndex);

				// Compute all rates and store them in your fields
				computeUptakeRate(s, mass, 0);
				//computeUptakeRate(s, mass);

				// Apply these rates on the grid
				fitUptakeRatesOnGrid(reacGrid, diffReacGrid, anAgent);
			}
		} else {
			for (ActiveAgent anAgent : _guild) {
				s = readConcentrationSeen(anAgent, concGrid);

				// Compute all rates and store them in your fields
				computeUptakeRate(s, anAgent);

				// Apply these rates on the grid
				fitUptakeRatesOnGrid(reacGrid, diffReacGrid, anAgent);
			}
		}

	}*/

	public abstract void updateMarginalMu(double[] s);
	public abstract double computeSpecRate(double[] s);
	public abstract double[] computeMarginalDiffMu(double[] s);

	/**
	 * @param reacGrid
	 * @param diffReacGrid
	 * @param anAgent
	 */
	public void fitUptakeRatesOnGrid(SoluteGrid[] reacGrid, SoluteGrid[] diffReacGrid,
			ActiveAgent anAgent) {
		DiscreteVector dC = new DiscreteVector(0, 0, 0);

		if (anAgent instanceof LocatedAgent) {

			for (int iGrid : _mySoluteIndex) {
				dC = reacGrid[iGrid].getDiscreteCoordinates(((LocatedAgent) anAgent).getLocation());
				reacGrid[iGrid].addValueAt(_uptakeRate[iGrid], dC);
				diffReacGrid[iGrid].addValueAt(_diffUptakeRate[iGrid], dC);
			}

		} else {
			for (int iGrid : _mySoluteIndex) {
				reacGrid[iGrid].addAllValues(_uptakeRate[iGrid]);
				diffReacGrid[iGrid].addAllValues(_diffUptakeRate[iGrid]);
			}
		}
	}

	/**
	 * Add the mass of all the agents of the guild on a received grid TODO : use
	 * this method for all ApplyReaction stuff ! The stored value is a
	 * CONCENTRATION
	 * @param aSpG
	 */
	public void fitAgentMassOnGrid(SpatialGrid aSpG) {
		for (ActiveAgent anActiveAgent : _guild) {
			anActiveAgent.fitMassOnGrid(aSpG, this._catalystIndex);
		}
	}

	/* ______________________ TOOLBOX _______________________________ */

	/**
	 * @param anAgent
	 * @return all the concentration seen by an agent on the default solute grid
	 */
	public double[] readConcentrationSeen(ActiveAgent anAgent, SoluteGrid[] concGrid) {

		double[] out = new double[concGrid.length];

		//sonia:chemostat
		//the concentration read by the agents is the one stored in the bulk (which has been previously updated)

		if (Simulator.isChemostat){



			for (int index =0; index<_soluteList.length; index++){


				for (AllBC aBC : _reacGrid.getDomain().getAllBoundaries()){
					if (aBC.hasBulk()){
						Bulk aBulk = aBC.getBulk();
						if(aBulk.getName().equals("chemostat")){
							out[index] = aBulk.getValue(_soluteList[index].soluteIndex);
						}
					}	
				}


				//System.out.println("solute " + _soluteList[index].getName());
				//System.out.println("concentrations seen by the agent...." + out[index]);

			}
		}else{

			if (anAgent instanceof LocatedAgent) {
				// The agent is a located agent, use the local concentration
				for (int iGrid = 0; iGrid<_soluteList.length; iGrid++){
					out[iGrid] = _soluteList[iGrid].getValueAround(((LocatedAgent) anAgent));
				}
			} else {
				// The agent is a non-located agent, use the average concentration
				for (int iGrid = 0; iGrid<_soluteList.length; iGrid++)
					out[iGrid] = _soluteList[iGrid].getAverage();
			}
		}

		return out;
	}

	public void writeReport(ResultFile bufferState, ResultFile bufferSum) throws Exception {
		fitGuildOnGrid();

		_reacGrid.writeReport(bufferState, bufferSum);
	}

	/**
	 * Build a grid with mass and mass-growth-rate of all agents of the guild
	 */
	public void fitGuildOnGrid() {

		_reacGrid.setAllValueAt(0);
		for (ActiveAgent anAgent : _guild) {
			// Sum mass of catalyser compartments on each grid cell
			anAgent.fitMassOnGrid(_guildGrid, _catalystIndex);

			// Apparent reaction rate on each grid cell
			anAgent.fitReacRateOnGrid(_reacGrid, reactionIndex);
		}
	}

	/* ______________ ACCESSORS _______________________________ */
	public double[] getSoluteYield() {
		return _soluteYield;
	}
	
	public int[] getSoluteIndex(){
		return _mySoluteIndex;
	}

	public double[] getParticulateYield() {
		return _particleYield;
	}


	public double[] getKinetic() {
		return _kineticParam;
	}

	public LinkedList<ActiveAgent> getGuild() {
		return _guild;
	}

	public abstract double[] getMarginalDiffMu();

}
