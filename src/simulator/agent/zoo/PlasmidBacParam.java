package simulator.agent.zoo;

import java.awt.Color;
import java.util.ArrayList;

import simulator.Simulator;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief TODO
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidBacParam extends BactEPSParam
{
	/**
	 * The maximum growth rate this PlasmidBac can achieve. Parameter for
	 * growth dependency of Plasmid scan rate.
	 */
	public double maxGrowthRate = 1.0;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public double lowTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public double highTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Parameter for collision frequency in chemostat simulations.
	 */
	public double collisionCoeff = 0.1;
	
	/**
	 * Whether or not to scale scan probabilities by distance from the host
	 * (only applies in biofilm simulations).
	 */
	public boolean scaleScanProb = false;
	
	/**
	 * Colours for POV-Ray output.
	 * 
	 * TODO define these for multiple hosted plasmids.
	 */
	public Color dColor, tColor, rColor;
	
	/**
	 * A list of Plasmid species that could be hosted by PlasmidBacs of this
	 * species. Useful in writing reports for output.
	 */
	public ArrayList<String> potentialPlasmids;
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidBacParam()
	{
		
	}
	
	/**
	 * 
	 */
	@Override
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		
		double tempDbl;
		Boolean tempBool;
		String tempCol;
		/*
		 * Get the growth dependency parameters.
		 */
		tempDbl = getSpeciesParameterDouble("maxGrowthRate",
											aSpeciesRoot, speciesDefaults);
		maxGrowthRate = Double.isFinite(tempDbl) ? tempDbl : maxGrowthRate;
		
		
		tempDbl = getSpeciesParameterDouble("lowTonusCutoff",
											aSpeciesRoot, speciesDefaults);
		lowTonusCutoff = Double.isFinite(tempDbl) ? tempDbl : lowTonusCutoff;
		
		
		tempDbl = getSpeciesParameterDouble("highTonusCutoff", aSpeciesRoot,
															speciesDefaults);
		highTonusCutoff = Double.isFinite(tempDbl)? tempDbl : highTonusCutoff;
		
		
		tempBool = getSpeciesParameterBool("scaleScanProb", aSpeciesRoot,
															speciesDefaults);
		scaleScanProb = (tempBool == null) ? scaleScanProb : tempBool;
		if ( scaleScanProb && Simulator.isChemostat )
		{
			LogFile.writeLogAlways("Cannot scale scan probabilities by"+
				"distance in the chemostat! Setting scaleScanProb to false");
			scaleScanProb = false;
		}
		
		/*
		 * Get the colours for POV-Ray output.
		 */
		tempCol = getSpeciesParameterString("donorColor",
											aSpeciesRoot, speciesDefaults);
		tempCol = ( tempCol == null ) ? "white" : tempCol;
		dColor = utils.UnitConverter.getColor(tempCol);
		
		
		tempCol = getSpeciesParameterString("transconjugantColor",
											aSpeciesRoot, speciesDefaults);
		tempCol = ( tempCol == null ) ? "white" : tempCol;
		tColor = utils.UnitConverter.getColor(tempCol);
		
		
		tempCol = getSpeciesParameterString("recipientColor",
											aSpeciesRoot, speciesDefaults);
		tempCol = ( tempCol == null ) ? "white" : tempCol;
		rColor = utils.UnitConverter.getColor(tempCol);
		
		potentialPlasmids = new ArrayList<String>();
	}
	
	/**
	 * \brief Add the given name to the list of Plasmid species that could be
	 * hosted by PlasmidBacs of this species.
	 * 
	 * <p>Checks that the given name is not already on the list.</p>
	 * 
	 * @param name Plasmid species name.
	 */
	public void addPotentialPlasmidName(String name)
	{
		if ( ! potentialPlasmids.contains(name) )
			potentialPlasmids.add(name);
	}
}
