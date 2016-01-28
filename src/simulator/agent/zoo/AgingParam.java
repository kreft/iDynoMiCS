package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SpatialGrid;
import utils.LogFile;
import utils.XMLParser;

public class AgingParam extends BacteriumParam
{
	/**
	 * Damage segregation parameter. Default alpha = 0 means symmetric
	 * division, whereas alpha = 1 would mean completely asymmetric. 
	 */
	Double alpha = 0.0;
	
	/**
	 * Stochastic deviation from alpha.
	 */
	Double alphaDev = 0.0;
	
	/**
	 * Damage repair investment parameter. Default beta = 0 means no
	 * investment in repair.
	 */
	Double beta  = 0.0;
	
	/**
	 * Is investment in repair optimal? Default is false, so repair takes the value
	 * for beta defined above
	 */
	boolean isOptimalRepair = false;
	
	/**
	 * Damage repair maximum rate parameter. Default value is 1.
	 */
	Double rMax  = 1.0;
	
	/**
	 * Repair yield parameter. Default value of 0.8 means 4 units of active
	 * biomass are recovered by repairing 5 units of damaged biomass.
	 */
	Double repY  = 0.8;
	
	/**
	 * Toxicity parameter. Default value is False, which means that damage is not toxic.
	 */
	boolean isToxic = false;
	
	/**
	 * Pint parameter. Default value is false, which means that division occurs 
	 * when total biomass reaches a certain parameter.
	 */
	boolean isPint = false;
	/**
	 * Linear growth paramater. Default value is false, which means that growth is exponential.
	 */
	boolean isLinear = false;
	
	public AgingParam()
	{
		super();
	}
		
	
	@Override
	public void init(Simulator aSim, XMLParser aSpeciesRoot, XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		Double value;
		/*
		 * Damage segregation parameter, alpha.
		 */
		value = getSpeciesParameterDouble("alpha", aSpeciesRoot, speciesDefaults);
		if ( Double.isFinite(value) )
			if ( value <= 1.0 && value >= 0.0 )
				alpha = value;
			else
				LogFile.writeLogAlways("alpha must be between 0 and 1!");
		/*
		 * Stochastic deviation from alpha.
		 */
		value = getSpeciesParameterDouble("alphaDev", aSpeciesRoot, speciesDefaults);
		if ( Double.isFinite(value) )
			if ( value <= 0.5 && value >= 0.0 )
				alphaDev = value;
			else
				LogFile.writeLogAlways("alphaDev must be between 0 and 0.5!");
		/*
		 * Damage repair maximum rate parameter.
		 */
		value = getSpeciesParameterDouble("beta", aSpeciesRoot, speciesDefaults);
		if ( Double.isFinite(value) )
			if ( value <= 1.0 && value >= 0.0 )
				beta = value;
			else
				LogFile.writeLogAlways("beta must be between 0 and 1!");
		/*
		 * Damage repair maximum rate parameter.
		 */
		value = getSpeciesParameterTime("rMax", aSpeciesRoot, speciesDefaults);
		rMax = Double.isFinite(value) ? value : rMax;
		/*
		 * Repair yield parameter.
		 */
		value = getSpeciesParameterDouble("repY", aSpeciesRoot, speciesDefaults);
		repY = Double.isFinite(value) ? value : repY;
		/*
		 * Whether damage is toxic
		 */
		Boolean boolTemp = getSpeciesParameterBool("isToxic",
				aSpeciesRoot, speciesDefaults);
		isToxic = (boolTemp == XMLParser.nullBool) ? isToxic : boolTemp;
		/*
		 * Whether division occurs at threshold value of total or active biomass
		 */
		boolTemp = getSpeciesParameterBool("isPint",
				aSpeciesRoot, speciesDefaults);
		isPint = (boolTemp == XMLParser.nullBool) ? isPint : boolTemp;
		/*
		 * Whether growth is linear
		 */
		boolTemp = getSpeciesParameterBool("isLinear",
				aSpeciesRoot, speciesDefaults);
		isLinear = (boolTemp == XMLParser.nullBool) ? isLinear : boolTemp;
		/*
		 * Whether beta is optimised
		 */
		boolTemp = getSpeciesParameterBool("isOptimalRepair",
				aSpeciesRoot, speciesDefaults);
		isOptimalRepair = (boolTemp == XMLParser.nullBool) ? isOptimalRepair : boolTemp;
	}
}
