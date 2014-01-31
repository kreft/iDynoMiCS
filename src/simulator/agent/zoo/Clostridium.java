package simulator.agent.zoo;

import Jama.Matrix;
import simulator.reaction.Reaction;
import simulator.Simulator;
import simulator.SoluteGrid;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief A bacterium with gene regulation, based on Clostridium acetobutylicum.
 * 
 * This species is based on the description of C. acetobutylicum given in:
 * Jabbari et al. (2013). The putative influence of the agr operon upon
 * survival mechanisms used by Clostridium acetobutylicum. Mathematical
 * Biosciences. 243: 223-239.
 * 
 * 
 * As of 
 * 
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 *
 */
public class Clostridium extends GeneRegBac implements Cloneable
{
	
	
	/**
	 * 
	 */
	protected Boolean _isSpore;
	
	
	protected SoluteGrid[] _soluteList;
	
	
	/**
	 * The index of the solute named "acid" in _soluteList. Simpler to store it
	 * here than to look it up every time!
	 */
	protected int _acidIndex;
	
	/**
	 * The acid concentration in this cell's immediate vicinity.
	 */
	protected Double _acidConc;
	
	
	
	
	public Clostridium()
	{
		super();
		_speciesParam = new ClostridiumParam();
		
		_isSpore = false;
		
		_numProtTypes = 1;
		
		_proteinNames = new String[_numProtTypes];
		_proteinLevels = ExtraMath.newDoubleArray(_numProtTypes);
		
		_proteinNames[0] = "Kin";
		_proteinLevels[0] = 1.0;
	}
	
	
	public Object clone() throws CloneNotSupportedException
	{
		Clostridium out = (Clostridium) super.clone();
		
		if (this._isSpore)
			LogFile.writeLogAlways("Warning: spore is being cloned!");
		out._isSpore = false;
		
		
		
		out._soluteList = this._soluteList;
		out._acidConc = 0.0;
		
		
		return out;
	}
	
	/**
	 * 
	 * @param aSim Simulator 
	 * @param aSpeciesRoot XMLParser 
	 */
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot)
	{
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		_soluteList = aSim.soluteList;
		_acidIndex = aSim.getSoluteIndex("acid");
		
	}
	
	/**
	 * 
	 */
	protected void internalStep()
	{
		
		updateGeneRegulation();
		
		grow();
		
		updateSize();
		
		if (willDivide())
			divide();
		
		if (willDie())
		{
			this.death = "tooSmall";
			die(true);
		}
		
		if (willSporulate())
			sporulate();
		
	}
	
	/**
	 * 
	 */
	public void updateGeneRegulation()
	{
		try
		{
			// Nothing to do if this cell has sporulated.
			if ( _isSpore )
				return;

			updateLocalSolutes();
			
			// Tell the solver to look for me, as I have the correct acid concn
			_regulationSolver.setReferenceAgent(this);

			Matrix protein = _regulationSolver.arrayToMatrix(_proteinLevels);
			protein = _regulationSolver.solve(protein,
					_agentGrid.AGENTTIMESTEP,
					getSpeciesParam().rtol, getSpeciesParam().hmax);
			_proteinLevels = _regulationSolver.matrixToArray(protein);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "Clostridium.updateGeneRegulation()");
		}
	}
	
	/**
	 * Finds the current solute concentrations for this cell and stores them in
	 * _soluteConc. 
	 */
	public void updateLocalSolutes()
	{
		try
		{
			_acidConc = _soluteList[_acidIndex].getValueAround(this);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "Clostridium.updateLocalSolutes()");
		}
	}
	
	/**
	 * 
	 * Note that no growth will happen if the cell is a spore, as no reactions
	 * are active.
	 */
	public void updateGrowthRates()
	{
		super.updateGrowthRates();
	}
	
	/**
	 * 
	 */
	public void grow()
	{
		updateGrowthRates();
		
		// We adjust the particle masses after calculating all the deltaParticle
		// values so that the reactions occur simultaneously.
		for (int i = 0; i < particleMass.length; i++)
			particleMass[i] += deltaParticle[i];
	}
	
	/**
	 * \brief Decide whether or not to sporulate based on current conditions.
	 * 
	 * TODO Make this a time-delayed process, i.e. not instant. See
	 * LocatedAgent.willDivide() for an idea of how to do this.
	 */
	public Boolean willSporulate()
	{
		return _proteinLevels[0] > 2.0;
	}
	
	/**
	 * \brief Cell sporulates instantly, switching off all reactions and
	 * halting growth.
	 */
	public void sporulate()
	{
		// Set the Boolean flag.
		_isSpore = true;
		
		// This cell is no longer growing.
		_netGrowthRate = 0.0;
		_netVolumeRate = 0.0;
		
		// All reactions are switched off.
		for (Reaction aReac : allReactions)
			switchOffreaction(aReac);
	}
	
	/**
	 * 
	 */
	public Matrix calc1stDeriv(Matrix levels)
	{
		Matrix dYdT = new Matrix(_numProtTypes, 1, 0.0);
		Double rate = 0.0;
		
		try
		{
		// Example: dK/dT = (c_K * H)/(H + 1) - alpha * K
		rate = getSpeciesParam().cK * _acidConc;
		rate /= _acidConc + 1;
		rate -= getSpeciesParam().alpha * levels.get(0, 0);
		dYdT.set(0, 0, rate);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "Clostridium.calc1stDeriv()");
			System.exit(-1);
		}
		
		return dYdT;
	}
	
	/**
	 * 
	 */
	public Matrix calcJacobian(Matrix levels)
	{
		Matrix dFdY = new Matrix(_numProtTypes, _numProtTypes, 0.0);
	
		try
		{
			// Example: d(dK/dT)/dK = - alpha
			dFdY.set(0, 0, - getSpeciesParam().alpha);
		}
		catch (Exception e)
		{
			LogFile.writeError(e, "Clostridium.calc1stDeriv()");
		}
		
		return dFdY;
	}
	
	
	/**
	 * \brief Specifies the header of the columns of output information for
	 * this agent.
	 * 
	 * Adds nothing new to the header from LocatedAgent. 
	 * Used in creation of results files.
	 * 
	 * @return	String specifying the header of each column of results
	 * associated with this agent.
	 */
	public StringBuffer sendHeader()
	{
		StringBuffer tempString = super.sendHeader();

		for (int i = 0; i < _numProtTypes; i++)
			tempString.append(","+_proteinNames[i]);
		
		tempString.append(",isSpore");
		
		return tempString;
	}
	
	/**
	 * \brief Creates an output string of information generated on this
	 * particular agent.
	 * 
	 * Used in creation of results files.
	 * Writes the data matching the header file.
	 * 
	 * @return	String containing results associated with this agent.
	 */
	public StringBuffer writeOutput()
	{
		StringBuffer tempString = super.writeOutput();
		
		for (int i = 0; i < _numProtTypes; i++)
			tempString.append(","+_proteinLevels[i]);
		
		tempString.append(","+_isSpore);
		
		return tempString;
	}
	
	public ClostridiumParam getSpeciesParam()
	{
		return (ClostridiumParam) _speciesParam;
	}
}
