/**
 * \package reaction
 * \brief Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS
 * 
 * Package of classes used to model stoichiometric and kinetic reactions in iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
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
/**
 * \brief Class containing methods to create and query reactions that occur in the simulation
 * 
 * Reactions in which solutes or biomass are produced or consumed mediate the dynamics of the entire iDynoMiCS simulation. The REACTION 
 * markup in the XML file specifies the reactions that occur. These may be defined via both stoichiometric and kinetic forms, and both of 
 * these representations are used in definining reactions in the protocol file.
 * 
 * @version 1.2
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Kieran Alden (k.j.alden@bham.ac.uk), Centre for Systems Biology, University of Birmingham, UK
 */
public abstract class Reaction implements Serializable 
{
	/**
	 * Name of this reaction. Set from the relevant tag in the protocol file
	 */
	public String                     reactionName;
	
	/**
	 * Index to this reaction in the simulation dictionary
	 */
	public int                        reactionIndex;
	
	/**
	 * Holds an integer that references the simulation particle dictionary - notes the particle that catalyses this reaction
	 */
	public int                        _catalystIndex;

	/**
	 * Agents hosting this process
	 */
	protected LinkedList<ActiveAgent> _guild = new LinkedList<ActiveAgent>();

	/**
	 * Local copy of the solute grids used in this simulation. Used for efficiency purposes 
	 */
	protected SoluteGrid[]            _soluteList;
	
	/**
	 * A list of solute indices: if that solute's yield in this reaction is non-zero then it is on this list.
	 * Hence, this is the list of solutes which are affected by this reaction. 
	 **/
	public int[]                   _mySoluteIndex;
	
	/** 
	 * Boolean noting whether or not the catalytic particle is affected by the reaction 
	*/
	public boolean					autocatalytic;
	
	/**
	 *	Solute grid that stores reaction rate for this reaction 
	 */
	protected SoluteGrid              _reacGrid;
	
	/**
	 * Grid concerning the biomass that catalyses this reaction
	 */
	protected SoluteGrid	_guildGrid;
	
	/**
	 * Buffer vector used to communicate uptake rate to the solver, for each solute
	 */
	public double[]                   totalUptake;

	/**
	 *  Temporary storage used during computation of reaction rates
	 */
	protected double                  _specRate;
	
	/**
	 * Buffer vector used to communicate diffision uptake rate to the solver, for each solute
	 */
	protected double[]                _diffUptakeRate;
	
	/**
	 * Buffer vector used to communicate uptake rate to the solver, for each solute
	 */
	protected double[]                _uptakeRate;
	
	/**
	 * Array to store the solute yield for this reaction.  If a solute is not concerned by the current reaction, the yield will be zero
	 */
	protected double[]                _soluteYield;
	
	/**
	 * Array to store the uptake or production of each solute, by this reaction, for the whole domain [KA August 2013]
	 */
	public double[]				  _soluteProductionOrUptake;
	
	/**
	 * Array to store the kinetic factors involved in this reaction
	 */
	protected double[]                _kineticParam;
	
	/** 
	 * Array to store the particle yield from this reaction. If a solute is not concerned by the current reaction, the yield will be zero
	 */
	protected double[]                _particleYield;

	/**
	 * Array to store the names of the particles yielded in this reaction (sonia: 21-05-09)
	 */
	protected String[]             _particleNameYield ;

	/**
	 * dilution rate from associated bulk (sonia 21-04-10)
	 */
	public double Dil;

	/**
	 * Temporary variable to hold number of solutes in this specified simulation case
	 */
	static int                        nSolute;
	
	/*************************************************************************************************************************
	 * CLASS METHODS 
	 ************************************************************************************************************************/


	/**
	 * \brief Initialises a reaction object, giving the object the properties stated in the protocol file
	 * 
	 * Reactions in which solutes or biomass are produced or consumed mediate the dynamics of iDynoMiCS. This constructor initialises 
	 * a particular reaction that is stated in the protocol file. Kinetic parameters are managed by the subdefined method in the 
	 * children class
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param aReactionRoot	The XML object containing the definition of one reaction in the protocol file
	 * @see Simulator.createReaction()
	 */
	public void init(Simulator aSim, XMLParser aReactionRoot) 
	{
		// Set the name of the reaction
		reactionName = aReactionRoot.getAttribute("name");

		// Get the number of solutes that exist in this simulation run
		nSolute = aSim.soluteList.length;
		
		// Get the number of declared particles that exist in this simulation run
		int nParticulate = aSim.particleDic.size();

		// Create a simple array of all solutes
		_soluteList = aSim.soluteList;
		
		// Initialise the reaction grid using the solute list grid array
		_reacGrid = new SoluteGrid(_soluteList[0]);
		
		// Set the reaction grid name
		_reacGrid.gridName = reactionName+"-rate";

		// Guild grid - biomass that catalyses the reaction 
		_guildGrid = new SoluteGrid(_soluteList[0]);
		_guildGrid.gridName = reactionName+"-pop";

		// Initialize buffer arrays
		_uptakeRate = new double[nSolute];
		_diffUptakeRate = new double[nSolute];
		totalUptake = new double[nSolute];

		// Initialize arrays storing reaction parameters
		_soluteYield = new double[nSolute];
		_soluteProductionOrUptake = new double[nSolute];
		_particleYield = new double[nParticulate];
		_particleNameYield = new String[nParticulate];

		
		// Set the boundary conditions if this is a chemostat condition
		for (AllBC aBC : _reacGrid.getDomain().getAllBoundaries())
		{
			if (aBC.hasBulk())
			{
				Bulk aBulk = aBC.getBulk();
				if(aBulk.getName().equals("chemostat"))
				{
					Dil = aBulk._D;
				}
			}	
		}

		// Extract the yields for solutes and particulates from the XML file
		fillParameters(aSim, aReactionRoot);
	}
	
	/**
	 * \brief Use the reaction class to fill the parameters fields of the agent. Populate only the yield vector, the kinetic vector will be filled by the subdefined method in the child classes
	 * 
	 * Use the reaction class to fill the parameters fields of the agent. Populate only the yield vector, the kinetic vector will be 
	 * filled by the subdefined method in the child classes
	 * 
	 * @param anAgent	The ActiveAgent which parameters are being populated
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param xmlRoot	The XML object containing the definition of one reaction in the protocol file
	 * @see Simulator.createReaction()
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
	 * \brief End of the initialization procedure. Store the name and simulation dictionary index for this reaction
	 * 
	 * End of the initialization procedure. Store the name and simulation dictionary index for this reaction
	 * 
	 * @param aReactionName	Name of this reaction
	 * @param aSimulator	Integer reference to the readcion simulation dictionary
	 * @see Simulator.addReaction()
	 */
	public void register(String aReactionName, Simulator aSimulator) {
		reactionName = aReactionName;
		reactionIndex = aSimulator.getReactionIndex(aReactionName);
	}

	/**
	 * \brief Utilised by a solver instance to know the solutes affected by this reaction
	 * 
	 * Utility method called by a solver instance to know the solutes affected by this reaction
	 * 
	 * @return a list of solutes affected by this reaction
	 * @see DiffusionSolver.addReaction()
	 */
	public LinkedList<String> declareSolutes() 
	{
		LinkedList<String> affectedSolutes = new LinkedList<String>();
		for (int index : _mySoluteIndex) 
		{
			affectedSolutes.add(_soluteList[index].getName());
		}

		return affectedSolutes;
	}

	/**
	 * \brief Reads the reaction information (name, solute yield, particle yield etc) from the XML file and initialises this reaction
	 * 
	 * This method takes one reaction from the XML file and initialises a Reaction object with the properties stated in the protocol 
	 * file. This includes reaction name, what the reaction is catalysed by, reaction kinetics, and solute and particle yields
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param xmlRoot	The XML object containing the definition of one reaction in the protocol file
	 */
	public void fillParameters(Simulator aSim, XMLParser xmlRoot) 
	{

		// Get the factor that is catalysing the reaction
		String catalystName = xmlRoot.getAttributeStr("catalyzedBy");
		_catalystIndex = aSim.particleDic.indexOf(catalystName);

		// Populate yield for solutes
		double yield;
		
		// Get the yield information from the reaction tag of the protocol file
		XMLParser parser = new XMLParser(xmlRoot.getChildElement("yield"));
		
		// Now determine how much of each declared solute yields from this reaction
		// KA 28/03/13 - Simplified things here - later there was a loop to count the number of non-zero yielding solutes
		// stored in jSolute. No real point looping again so I've added that count in here. Then there was a further loop to cycle 
		// through this again, making an array of non-zero solutes. So I've scrapped this - set an array equal to the number of solutes,
		// and resized this after population such that it is the size of the number of non-zero yielding solutes.
		int nonZeroSolutes = 0;
		int[] _mySoluteIndexMax = new int[_soluteList.length];
		
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++) 
		{
			// Get the yield of this solute if present in the protocol file
			yield = parser.getParamDbl( _soluteList[iSolute].getName());
			if (!Double.isNaN(yield)) 
			{
				// Add the yield of this solute to the array
				_soluteYield[iSolute] = yield;
				
				// KA - Add to the array of solutes that are non-zero yielding
				_mySoluteIndexMax[nonZeroSolutes] = iSolute;
						
				// KA - Increase the count of the solutes with a non-zero yield
				nonZeroSolutes++;
				
			} 
			else 
			{
				// If not present assume a yield of 0
				_soluteYield[iSolute] = 0.0;
			}
		}
		
		// Now resize the non-zero solutes array
		_mySoluteIndex = resizeArray(_mySoluteIndexMax,nonZeroSolutes);
		
		
		
		// Populate yields for particles from this reaction
		String particleName;
		
		// Go through  each of the particles declared in this simulation and establish if any yield from reaction
		for (int iParticle = 0; iParticle<aSim.particleDic.size(); iParticle++) 
		{
			// Get the yield from the XML file, if present
			yield = parser.getParamDbl(aSim.particleDic.get(iParticle));
			
			// Get the name of the particle from the dictionary
			particleName = aSim.particleDic.get(iParticle);
			if (!Double.isNaN(yield)) 
			{
				// If there is a yield, store in the particle yield array
				_particleYield[iParticle] = yield;
				
				// Store the name of the particle that is yielded in the second array
				_particleNameYield[iParticle] = particleName;
			}
		}
		
		// Determine whether the catalytic particle is affected by the reaction,
		// and hence whether the reaction is autocatalytic or not.
		// This is important in ActiveAgent.grow()
		
		if (_particleYield[_catalystIndex]==0)
		{
			autocatalytic = false;
			LogFile.writeLog("Reaction "+reactionIndex+" is not autocatalyic");
		} 
		else 
		{
			autocatalytic = true;
			LogFile.writeLog("Reaction "+reactionIndex+" is autocatalyic");
		}
	}

	/* ______________ METHODS FOR REACTION MANAGEMENT_________________________ */

	/**
	 * \brief Register an agent among the guild of this pathway
	 *
	 * Register an agent among the guild of this pathway
	 * 
	 * @param anAgent	ActiveAgent to register
	 */
	public void addAgent(ActiveAgent anAgent) {
		_guild.add(anAgent);
	}

	/**
	 * \brief Remove an agent among the guild of this pathway
	 * 
	 * Remove an agent among the guild of this pathway
	 * 
	 * @param anAgent	ActiveAgent to remove
	 */
	public void removeAgent(ActiveAgent anAgent) {
		_guild.remove(anAgent);
	}

	/**
	 * \brief Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 *  Compute the marginal growth rate (i.e the specific growth rate times the mass of the particle which is mediating this reaction)
	 * 
	 * @param anAgent	Specific growth rate for this ActiveAgent
	 * @return	The marginal growth rate
	 */
	public abstract double computeMassGrowthRate(ActiveAgent anAgent);
	
	/**
	 * \brief Compute the specific growth rate
	 * 
	 * Compute the specific growth rate
	 * 
	 * @param anAgent	Specific growth rate for this ActiveAgent
	 * @return	The specific growth rate
	 */
	public abstract double computeSpecGrowthRate(ActiveAgent anAgent);

	/**
	 * \brief Return the specific reaction rate for a given agent
	 * 
	 * Return the specific reaction rate for a given agent
	 * 
	 * @param anAgent	Agent to use to determine solute concentration and calculate reaction rate
	 * @see ActiveAgent.grow()
	 * @see Episome.computeRate(EpiBac)
	 */
	public abstract void computeSpecificGrowthRate(ActiveAgent anAgent);

	/**
	 * \brief Compute specific growth rate in function to concentrations sent
	 * 
	 * Compute specific growth rate in function to concentrations sent
	 * 
	 * @param s	Array of solute concentration
	 * @param anAgent	Parameters used are those defined for THIS agent
	 */
	public abstract void computeSpecificGrowthRate(double[] s, ActiveAgent anAgent);

	/**
	 * \brief Compute specific growth rate in function of concentrations sent Parameters used are those defined for the reaction.
	 * 
	 * Compute specific growth rate in function of concentrations sent Parameters used are those defined for the reaction.
	 * 
	 * @param s	Array of solute concentration
	 */
	public abstract void computeSpecificGrowthRate(double[] s);

	/**
	 * \brief Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * Update the array of uptake rates and the array of its derivative. Based on default values of parameters. Unit is fg.h-1
	 * 
	 * @param s	The concentration of solute locally observed
	 * @param mass	Mass of the catalyst (cell...)
	 * @param t	Time
	 */
	public abstract void computeUptakeRate(double[] s, double mass, double t);

	/**
	 * \brief Compute reaction rate on each concerned solute grids Assumes same parameters for all the agents of a same guild
	 * 
	 * Compute reaction rate on each concerned solute grids Assumes same parameters for all the agents of a same guild
	 * 
	 * @param concGrid	Solute concentration
	 * @param reacGrid	Contribution of the reaction to the solute concentration dynamics
	 * @param diffReacGrid	Derivative of the previous grid
	 * @param biomassGrid	CONCENTRATION of the catalyst
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
					computeUptakeRate(s, biomassGrid.grid[i][j][k],0);

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

	

	/**
	 * \brief Add the mass of all the agents of the guild on a received grid. The stored value is a CONCENTRATION
	 * 
	 * Add the mass of all the agents of the guild on a received grid. The stored value is a CONCENTRATION
	 * 
	 * @param aSpG	Spatial grid of reactive biomass
	 */
	public void fitAgentMassOnGrid(SpatialGrid aSpG) {
		for (ActiveAgent anActiveAgent : _guild) {
			anActiveAgent.fitMassOnGrid(aSpG, this._catalystIndex);
		}
	}

	/**
	 * \brief Return a double array with all concentration seen by an agent on the default solute grid
	 * 
	 * Return a double array with all concentration seen by an agent on the default solute grid
	 * 
	 * @param anAgent	ActiveAgent responding to solute levels
	 * @param concGrid	Solute concentration grid
	 * @return all the concentration seen by an agent on the default solute grid
	 */
	public double[] readConcentrationSeen(ActiveAgent anAgent, SoluteGrid[] concGrid) {

		double[] out = new double[concGrid.length];

		//sonia:chemostat
		//the concentration read by the agents is the one stored in the bulk (which has been previously updated)
		if (Simulator.isChemostat)
		{
			for (int index =0; index<_soluteList.length; index++)
			{
				for (AllBC aBC : _reacGrid.getDomain().getAllBoundaries())
				{
					if (aBC.hasBulk())
					{
						Bulk aBulk = aBC.getBulk();
						if(aBulk.getName().equals("chemostat"))
						{
							out[index] = aBulk.getValue(_soluteList[index].soluteIndex);
						}
					}	
				}
			}
		}
		else
		{
			if (anAgent instanceof LocatedAgent) 
			{
				// The agent is a located agent, use the local concentration
				for (int iGrid = 0; iGrid<_soluteList.length; iGrid++)
				{
					out[iGrid] = _soluteList[iGrid].getValueAround(((LocatedAgent) anAgent));
				}
			} 
			else 
			{
				// The agent is a non-located agent, use the average concentration
				for (int iGrid = 0; iGrid<_soluteList.length; iGrid++)
					out[iGrid] = _soluteList[iGrid].getAverage();
			}
		}

		return out;
	}
	
	/**
	 * \brief Calculates the production or uptake of each solute for this reaction, over the whole domain, for this simulation timestep
	 * 
	 * Calculates the production or uptake of each solute for this reaction, over the whole domain, for this simulation timestep. These 
	 * are totalled for all reactions and written to the env_state files at the appropriate output period
	 */
	public void calculateSoluteChange()
	{
		// For each solute in the simulation
		for (int iSolute = 0; iSolute<_soluteList.length; iSolute++) 
		{
			// To obtain figure for this solute multiply the sum of the reaction rates by the stochiometry
			_soluteProductionOrUptake[iSolute] = _soluteYield[iSolute]*getReactionRateSum();
			
		}
	}

	/**
	 * \brief Write reaction information to the result file stream
	 * 
	 * Write reaction information to the result file stream
	 * 
	 * @param bufferState	Buffer for state result files
	 * @param bufferSum	Buffer for sum result files
	 * @throws Exception	Exception thrown if these files cannot be opened
	 */
	public void writeReport(ResultFile bufferState, ResultFile bufferSum) throws Exception {
		fitGuildOnGrid();

		_reacGrid.writeReport(bufferState, bufferSum);
	}

	/**
	 * \brief Build a grid with mass and mass-growth-rate of all agents of the guild
	 * 
	 * Build a grid with mass and mass-growth-rate of all agents of the guild. Used when outputting simulation state to file
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
	
	/**
	* \brief Utility to resize an array. Used in case of creating an array of non-zero solutes
	* 
	* The former iDynoMiCS included three different loops that went through all solutes to find those that are non-zero yielding. This 
	* seemed pointless and has been reduced to one. However one of these builds an array of solutes that are non-zero yielding and thus 
	* needed an answer from the first loop. Thus, the array is now initialised at the max size it could be (number of all solutes) and 
	* populated when a non-zero solute is found. This function is then used to resize the array and get rid of any non-necessary space
	* 
	* @added 28/03/13 - Kieran Alden
	* @param oldArray  the old array, to be reallocated.
	* @param newSize   the new array size.
	* @return          A new array with the same contents.
	*/
	private static int[] resizeArray (Object oldArray, int newSize) 
	{
	   int oldSize = java.lang.reflect.Array.getLength(oldArray);
	   int[] newArray = new int[newSize];
	   int preserveLength = Math.min(oldSize, newSize);
	   if (preserveLength > 0)
	      System.arraycopy(oldArray, 0, newArray, 0, preserveLength);
	   return newArray; 
	}

	/**
	 * \brief Return the solute yield
	 * 
	 * Return the solute yield
	 * 
	 * @return Double array of solute yield
	 */
	public double[] getSoluteYield() {
		return _soluteYield;
	}
	
	/**
	 * \brief Return the solute index
	 * 
	 * Return the solute index
	 * 
	 * @return Integer array of solute index
	 */
	public int[] getSoluteIndex(){
		return _mySoluteIndex;
	}

	/**
	 * \brief Return the particulate yield
	 * 
	 * Return the particulate yield
	 * 
	 * @return Double array of particulate yield
	 */
	public double[] getParticulateYield() {
		return _particleYield;
	}


	/**
	 * \brief Return the kinetic parameters used in this reaction
	 * 
	 * Return the kinetic parameters used in this reaction
	 * 
	 * @return	Double array containing kinetic parameters
	 */
	public double[] getKinetic() {
		return _kineticParam;
	}

	/**
	 * \brief Return the guild of this pathway
	 * 
	 * Return the guild of this pathway
	 * 
	 * @return	LinkedList containing all agents in the guild
	 */
	public LinkedList<ActiveAgent> getGuild() {
		return _guild;
	}
	
	/**
	 * \brief Calculate the rate of change of each uptake rate with respect to each solute. Returned as a matrix
	 * 
	 * Calculate the rate of change of each uptake rate with respect to each solute. Returned as a matrix
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param biomass	Total particle mass in the system which catalyses this reaction
	 * @return Matrix containing rate of change of each uptake rate with respect to each solute
	 */ 
	public abstract Matrix calcdMUdS(Matrix S, double biomass);
	
	/**
	 * \brief Calculate the rate of change of each uptake rate with respect to time.
	 * 
	 * Calculate the rate of change of each uptake rate with respect to time. dMUdT = catalyticBiomass*specificGrowthRate*soluteYield.
	 * Returned as a matrix
	 * 
	 * @param S	Temporary container for solute concentration
	 * @param biomass	Total particle mass in the system which catalyses this reaction
	 * @return Matrix containing rate of change of each uptake rate with respect to time
	 */ 
	public abstract Matrix calcdMUdT(Matrix S, double biomass);

	/**
	 * \brief Update the Marginal Mu data matrix
	 * 
	 * Update the Marginal Mu data matrix
	 * 
	 * @param s	Temporary container for solute concentration 
	 */
	public abstract void updateMarginalMu(double[] s);
	
	/**
	 * \brief Compute the specific growth rate
	 * 
	 * Compute the specific growth rate. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return	The specific growth rate
	 */
	public abstract double computeSpecRate(double[] s);
	
	/**
	 * \brief Compute the marginal difference array
	 * 
	 * Compute the marginal difference array. Don't forget to update marginalMu before calling this! 
	 * 
	 * @param s	Temporary container for solute concentration 
	 * @return Marginal diff array
	 */
	public abstract double[] computeMarginalDiffMu(double[] s);
	
	public double getReactionRateSum()
	{
		return _reacGrid.getSum();
	}


}
