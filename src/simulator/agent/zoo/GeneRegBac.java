package simulator.agent.zoo;

import Jama.Matrix;
import odeSolver.GeneRegSolver;
import simulator.Simulator;
import utils.XMLParser;

public abstract class GeneRegBac extends Bacterium implements Cloneable
{
	/**
	 * Ordinary Differential Solver used to calculate the changes in gene
	 * expression.
	 */
	protected GeneRegSolver _regulationSolver = new GeneRegSolver();
	
	/**
	 * 
	 */
	protected int _numProtTypes;
	
	/**
	 * 
	 */
	protected String[] _proteinNames;
	
	/**
	 * 
	 */
	protected Double[] _proteinLevels;
	
	
	
	public GeneRegBac()
	{
		super();
		_speciesParam = new GeneRegBacParam();
	}
	
	
	public Object clone() throws CloneNotSupportedException
	{
		GeneRegBac out = (GeneRegBac) super.clone();
		out._regulationSolver = this._regulationSolver;
		out._numProtTypes = this._numProtTypes;
		out._proteinNames = this._proteinNames;
		out._proteinLevels = this._proteinLevels;
		return out;
	}
	
	
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot)
	{
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		_regulationSolver.init(_numProtTypes);
		_regulationSolver.setReferenceAgent(this);
	}
	
	
	public abstract Matrix calc1stDeriv(Matrix levels);
	
	public abstract Matrix calcJacobian(Matrix levels);
	
	
	public GeneRegBacParam getSpeciesParam()
	{
		return (GeneRegBacParam) _speciesParam;
	}
	
}