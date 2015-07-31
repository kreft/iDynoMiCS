package simulator.agent.zoo;

import java.awt.Color;
import java.util.ArrayList;

import simulator.Simulator;
import utils.XMLParser;

/**
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidBacParam extends BactEPSParam
{
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public Double lowTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public Double highTonusCutoff = - Double.MAX_VALUE;
	
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
		
		Double tempDbl;
		String tempCol;
		/*
		 * Get the growth dependency parameters.
		 */
		tempDbl = getSpeciesParameterDouble("lowTonusCutoff",
											aSpeciesRoot, speciesDefaults);
		lowTonusCutoff = Double.isFinite(tempDbl) ? tempDbl : lowTonusCutoff;
		
		tempDbl = getSpeciesParameterDouble("highTonusCutoff", aSpeciesRoot,
															speciesDefaults);
		highTonusCutoff = Double.isFinite(tempDbl) ? tempDbl : highTonusCutoff;
		
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
