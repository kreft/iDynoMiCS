/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * top-level class of the simulation core. It is used to create and run a
 * simulation; this class is called by the class Idynomics.
 * 
 */

/**
 * 
 * @since June 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology, University of Birmingham (UK)
 * 
 */

package simulator.agent;

import java.util.List;
import java.io.Serializable;
import java.awt.Color;

import org.jdom.Element;

import simulator.Simulator;
import simulator.geometry.*;
import utils.LogFile;
import utils.XMLParser;
import utils.ExtraMath;

public class Species implements Serializable {

	// Serial version used for the serialisation of the class
	private static final long  serialVersionUID = 1L;

	public Simulator           currentSimulator;
	public String              speciesName;
	public int                 speciesIndex;
	public Color               color;

	// Coomputation domain of the species
	public Domain              domain;

	protected SpecialisedAgent _progenitor;
	protected int              _population      = 0;

	/* ___________________ CONSTRUCTOR _____________________________________ */
	public Species(Simulator aSimulator, XMLParser aSpRoot) {
		speciesName = aSpRoot.getAttribute("name");
		String colorName = aSpRoot.getParam("color");
		if (colorName==null) colorName = "white";
		color = utils.UnitConverter.getColor(colorName);

		// Register it
		speciesIndex = aSimulator.speciesList.size();
		currentSimulator = aSimulator;

		// Create the progenitor and tune its speciesParam object
		_progenitor = (SpecialisedAgent) aSpRoot.instanceCreator("simulator.agent.zoo");
		_progenitor.getSpeciesParam().init(aSimulator, aSpRoot);
		_progenitor.setSpecies(this);

		domain = aSimulator.world.getDomain(aSpRoot.getParam("computationDomain"));
	}

	public Species(SpecialisedAgent aProgenitor) {
		_progenitor = aProgenitor;
		aProgenitor.setSpecies(this);
	}

	/**
	 * Register the created species to speciesManager
	 * @param aSpecies
	 */
	public void register(Simulator aSimulator, XMLParser aSpRoot) {
		currentSimulator = aSimulator;
		speciesName = aSpRoot.getAttribute("name");

		speciesIndex = aSimulator.speciesList.size();
		aSimulator.speciesList.add(this);

		domain = aSimulator.world.getDomain(aSpRoot.getAttribute("computationDomain"));
	}

	/**
	 * Create clones of the progenitor within the birth area specified before
	 * @param XML markup
	 */
	public void createPop(XMLParser spRoot) {

		double howMany = spRoot.getAttributeDbl("number");
		if (howMany==0) return;

		// Define the birth area
		ContinuousVector[] _initArea = defineSquareArea(spRoot);

		// Create all the required agents
		ContinuousVector cc;
		for (int i = 0; i<howMany; i++) {
			if (_progenitor instanceof LocatedAgent) {
				cc = new ContinuousVector();
				//sonia:chemostat
				if(Simulator.isChemostat){
					//do not randomly generate coordinates for the agent within the iniArea given in the protocol file
				}else{
					shuffleCoordinates(cc, _initArea);
				}
				((LocatedAgent) _progenitor).createNewAgent(cc);
			} else {
				_progenitor.createNewAgent();
			}
		}

		LogFile.writeLog(howMany+" agents of species "+speciesName+" successfully created");
	}

	public void notifyBirth() {
		_population++;
	}

	public void notifyDeath() {
		_population--;
	}

	/**
	 * @return a clone of the progenitor
	 * @throws CloneNotSupportedException
	 */
	public SpecialisedAgent sendNewAgent() throws CloneNotSupportedException {
		return _progenitor.sendNewAgent();
	}

	/* ______________ TOOLS ____________________________________________ */

	public int getPopulation() {
		return _population;
	}

	public SpecialisedAgent getProgenitor() {
		return _progenitor;
	}

	public SpeciesParam getSpeciesParam() {
		return _progenitor.getSpeciesParam();
	}

	public Species getSpecies(String speciesName) {
		return currentSimulator.speciesList.get(currentSimulator.getSpeciesIndex(speciesName));
	}

	public ContinuousVector[] defineSquareArea(XMLParser spRoot) {
		List<Element> area = spRoot.getChildren("coordinates");
		ContinuousVector cc1 = new ContinuousVector((Element) area.get(0));
		ContinuousVector cc2 = new ContinuousVector((Element) area.get(1));

		ContinuousVector[] initArea = new ContinuousVector[2];
		initArea[0] = new ContinuousVector();
		initArea[1] = new ContinuousVector();

		initArea[0].x = Math.min(cc1.x, cc2.x);
		initArea[0].y = Math.min(cc1.y, cc2.y);
		initArea[0].z = Math.min(cc1.z, cc2.z);
		initArea[1].x = Math.max(cc1.x, cc2.x);
		initArea[1].y = Math.max(cc1.y, cc2.y);
		initArea[1].z = Math.max(cc1.z, cc2.z);

		// Inthe case of 2D simulation, the agent's z-coordinate is 0.
		if (!domain.is3D) {
			initArea[0].z = 0;
			initArea[1].z = 0;
		}
		return initArea;
	}

	/**
	 * Generate random but valid continuous coordinates inside a volume
	 * @param cc
	 * @param area
	 */
	public void shuffleCoordinates(ContinuousVector cc, ContinuousVector[] area) {
		boolean test = true;
		while (test) {
			cc.x = area[0].x+ExtraMath.getUniRand()*(area[1].x-area[0].x);
			cc.y = area[0].y+ExtraMath.getUniRand()*(area[1].y-area[0].y);
			cc.z = area[0].z+ExtraMath.getUniRand()*(area[1].z-area[0].z);
			test = !(domain.testCrossedBoundary(cc)==null);

		}
	}
}
