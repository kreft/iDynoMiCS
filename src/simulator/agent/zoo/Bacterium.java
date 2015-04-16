/**
 * \package simulator.agent.zoo
 * \brief Package of agents that can be included in iDynoMiCS and classes to
 * store parameters for these agent types.
 * 
 * Package of agents that can be included in iDynoMiCS and classes to store
 * parameters for these agent types. This package is part of iDynoMiCS v1.2,
 * governed by the CeCILL license under French law and abides by the rules of
 * distribution of free software. You can use, modify and/ or redistribute
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS
 * and INRIA at the following URL  "http://www.cecill.info".
 */
package simulator.agent.zoo;

import java.awt.Color;
import java.math.BigInteger;

import simulator.agent.*;
import simulator.geometry.ContinuousVector;
import simulator.Simulator;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Creates an object of the Bacterium Species. Acts as a superclass for
 * all current species with the exception of ParticulateEPS.
 * 
 * The first part of the Bacterium species mark-up in the XML file lists the
 * different particles composing this species; these particles form
 * 'compartments' within the agent that signify the types of biomass that the
 * agent contains. The possible particle types that may be used are taken from
 * the particle mark-ups defined early in the protocol file, generally
 * consisting of 'biomass', 'inert', and 'capsule'. The initial mass of each
 * particle may be given (in units of femtograms, 1 fg = 10-15 g), but if you
 * set the value to zero the simulator will assign a reasonable random initial
 * value for the mass. For this Bacterium species, only one particle type is 
 * actually necessary to include, but it is common to include several types. In
 * this species, inert and capsule compartments exhibit special behaviour: 
 * inert biomass will accumulate in the agent, while a capsule will be excreted
 * if it accumulates to take up too much of the agent's mass. Note that the
 * 'capsule' particle also requires you to specify its type via the class
 * attribute; when the capsule has become too large and needs to be excreted,
 * an agent of that class will be created (this means that it is necessary to 
 * have previously defined a species with that name.
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany).
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France.
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of
 * Engineering Sciences and Applied Mathematics, Northwestern University (USA).
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 *
 */
public class Bacterium extends LocatedAgent implements Cloneable 
{	
	/**
	 * Boolean noting whether any EPS particles have been declared as part of
	 * this Bacterium.
	 */
	protected Boolean _hasEps = false;
	
	/**
	 * Boolean noting whether any inert particles have been declared as part
	 * of this Bacterium.
	 */
	protected Boolean _hasInert = false;
	
	/**
	 * Previously declared species object of the type of agent that will be
	 * created when the capsule is excreted.
	 */
	protected Species _epsSpecies;

	/**
	 * \brief Constructor used to generate progenitor and initialise an object
	 * to store relevant parameters
	 */
	public Bacterium() 
	{
		super();
		_speciesParam = new BacteriumParam();
	}

	/**
	 * \brief Used by makeKid method to create a daughter Bacterium agent by
	 * cloning this agent and parameter objects.
	 * 
	 * @throws CloneNotSupportedException 	Thrown if the agent cannot be
	 * cloned.
	 */
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		Bacterium out = (Bacterium) super.clone();
		out._hasEps = this._hasEps;
		out._epsSpecies = this._epsSpecies;
		return out;
	}

	/**
	 * \brief Creates a Bacterium agent from the parameters specified in the
	 * XML protocol file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aSpeciesRoot	A species mark-up within the specified protocol
	 * file.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot) 
	{
		// Initialisation of the Active agent
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		// Species index 
		int spIndex;
		
		// Check if it is a EPS-producing species
		for (XMLParser parser : aSpeciesRoot.getChildrenParsers("particle"))
		{
			if (parser.getName().equals("capsule"))
			{
				_hasEps = true;
				spIndex = aSim.getSpeciesIndex(parser.getAttribute("class"));
				_epsSpecies = aSim.speciesList.get(spIndex);
			}
			if (parser.getAttribute("name").equals("inert"))
				_hasInert = true;			
		}

		/* If no mass defined, use the division radius to find the mass */
		if ( this._totalMass.equals(0.0))
		{
			guessMass();
			LogFile.writeLog("Guessing "+this.getSpecies().speciesName+
										" initial mass at: "+this._totalMass);
		}

		// SET CELL RADIUS, GENERATION, AND GENEOLOGY
		init();
	}

	/**
	 * \brief Create an agent using information in a previous state or
	 * ]initialisation file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file
	 * that is used to recreate this agent.
	 */
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
	{
		// this writes no unique values, so doesn't need unique reading-in
		// (for a template on how to read in data, look in LocatedAgent.java)
		super.initFromResultFile(aSim,singleAgentData);
	}

	/**
	 * \brief Initialises any new agent (progenitor or daughter cell), setting
	 * cell radius, generation, and genealogy.
	 */
	public void init() 
	{
		// Lineage management : this is a new agent, he has no known parents.
		_generation = 0;
		_genealogy = BigInteger.ZERO;
		// Determine the radius, volume and total mass of the agent.
		updateSize();
	}

	/**
	 * \brief Called by Bacterium.createAgent and to obtain another instance
	 * of the same species (totally independent).
	 * 
	 * The returned agent is NOT registered.
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object
	 * cannot be cloned
	 */
	@Override
	public Bacterium sendNewAgent() throws CloneNotSupportedException 
	{
		// Clone the agent and initialise it
		Bacterium baby = (Bacterium) this.clone();
		baby.init();
		return baby;
	}

	/**
	 * \brief Create a new Bacterium agent (who a priori is registered in at
	 * least one container).
	 * 
	 * This agent is located on the relevant grid.
	 */
	@Override
	public void createNewAgent(ContinuousVector position) 
	{
		try 
		{
			// Get a clone of the progenitor
			Bacterium baby = sendNewAgent();
			baby.giveName();

			updateMass();
			/* If no mass defined, use the division radius to find the mass */
			// Note this should have been done already in initFromProtocolFile
			if ( this._totalMass == 0.0 )
			{
				guessMass();
				LogFile.writeLog("Warning: Bacterium.createNewAgent calling guessMass()");
			}
			
			// randomise its mass
			baby.randomiseMass();
			//System.out.println("RADIUS AT THIS POINT: "+this._totalRadius);
			baby.updateSize();
			//System.out.println("RADIUS AFTER CALL: "+this._totalRadius);

			// Just to avoid to be in the carrier
			position.x += this._totalRadius;
			
			baby.setLocation(position);
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			utils.LogFile.writeLog("Error met in Bacterium.createNewAgent(): "+e);
		}
	}

	/**
	 * \brief Mutates inherited parameters and distributes particle mass -
	 * either exponentially or normally, dependent on value of distMethod.
	 */
	public void randomiseMass() 
	{	
		// distMethod true -> distribute exponentially
		// distMethod false -> distribute normally
		if (getSpeciesParam().distMethod)
			for (int i = 0; i < particleMass.length; i++)
				particleMass[i] *= ExtraMath.getExp2Rand();
		else
			for (int i = 0; i<particleMass.length; i++)
				particleMass[i] = ExtraMath.deviateFromCV(
						1.5*particleMass[i], getSpeciesParam().initialMassCV);
	}
	
	/**
	 * \brief Used by Bacterium.divide() method to create a daughter cell of
	 * this agent.
	 * 
	 * @throws CloneNotSupportedException	Exception thrown if the object
	 * cannot be cloned.
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		super.makeKid();
	}

	/* ___________________ STEP METHODS _______________________________ */
	
	/**
	 * \brief Called at each time step of the simulation to compute mass
	 * growth and update radius, mass, and volume.
	 * 
	 * Under the control of the method Step of the class Agent to avoid
	 * multiple calls. Also determines whether the agent has reached the size
	 * at which it must divide, and monitors agent death.
	 */
	@Override
	protected void internalStep()
	{
		/*
		 * Compute mass growth over all compartments.
		 */
		grow();
		updateSize();
		/*
		 * Test if the EPS capsule has to be excreted.
		 */
		manageEPS();
		/*
		 * Divide if you have to.
		 */
		if ( willDivide() )
			divide();
		/*
		 * Die if you have to.
		 */
		if ( willDie() )
		{
			this.death = "tooSmall";
			die(true);
		}
	}
	
	/**
	 * \brief Converts the agent into particle EPS and inert on agent death.
	 * 
	 * @param isStarving	Boolean noting whether the agent currently has
	 * access to any resources.
	 */
	@Override
	public void die(Boolean isStarving)
	{
		super.die(isStarving);
		if ( isStarving && _hasEps )
			excreteEPS(1.0);
	}

	/**
	 * \brief Excretes EPS particle with 75% of initial mass if the EPS
	 * capsule is too thick.
	 * 
	 * TODO Make it possible to set these parameters in the protocol file.
	 */
	public void manageEPS()
	{
		if ( ! _hasEps )
			return;
		/*
		 * Manage excretion.
		 */
		if ( _volume < _totalVolume*(1-getSpeciesParam().epsMax) )
			excreteEPS(ExtraMath.getUniRand(.6, .9));
	}

	/**
	 * \brief Excretes part of the agents bound EPS as slime EPS.
	 * 
	 * @param ratio	Ratio of EPS that should be excreted.
	 */
	public void excreteEPS(Double ratio)
	{
		int indexEPS = getIndexCapsule();
		/*
		 * Check that the mass to excrete does exist.
		 */
		if ( particleMass[indexEPS]*ratio <= 0.0 )
			return;
		/*
		 * Create the particle.
		 */
		ParticulateEPS eps = (ParticulateEPS) _epsSpecies.getProgenitor();
		/*
		 * If the particle has been successfully created, update your size.
		 */
		if ( eps.createByExcretion(this, ratio) )
		{
			particleMass[indexEPS] *= (1.0 - ratio);
			updateSize();
		}
	}
	
	/**
	 * \brief Sets the mass of the primary particle of the progenitor to half
	 * the mass at-division (i.e. the average mass after division, regardless
	 * of babyMassFrac).
	 * 
	 * The mass-at-division is determined by the particle density and the
	 * volume-at-division; volume is determined by the divRadius, or the
	 * divRadius and the z-resolution if in a 2D simulation.
	 */
	public void guessMass()
	{
		Double val = getSpeciesParam().divRadius;
		/*
		 * We calculate the mass-at-division:
		 *  - in chemostats and 3D the cell is spherical.
		 *  - in 2D the cell is cylindrical.
		 */
		if (Simulator.isChemostat || _agentGrid.is3D)
			val = ExtraMath.volumeOfASphere(val); 
		else
			val = ExtraMath.volumeOfACylinder(val, _species.domain.length_Z);
		particleMass[0] = getSpeciesParam().particleDensity[0] * val / 2;
		updateMass();
	}
	
	/* _______________ RADIUS, MASS AND VOLUME _____________________ */

	/**
	 * \brief Update the volume of this agent by examining the particle
	 * density.
	 */
	@Override
	public void updateVolume()
	{
		Double[] density = getSpeciesParam().particleDensity;
		_totalVolume = 0.0;
		for ( int i = 0; i < particleMass.length; i++ )
			_totalVolume += particleMass[i] / density[i];
  		/*
  		 * Compute volume with or without EPS capsule.
  		 */
		_volume = _totalVolume;
		if ( _hasEps )
		{
			int i = getIndexCapsule();
  			_volume -= (particleMass[i] / density[i]);
  		}
  	}
	
	/**
	 * \brief Return the set of parameters that is associated with the object
	 * of this species.
	 * 
	 * @return Object of BacteriumParam that stores the parameters associated
	 * with this species.
	 */
	@Override
	public BacteriumParam getSpeciesParam()
	{
		return (BacteriumParam) _speciesParam;
	}

	/**
	 * \brief Return whether this bacterium contains any EPS particles
	 * (capsule in the protocol file).
	 * 
	 *  @return Boolean noting whether this bacterium object contains EPS
	 *  particles.
	 */
	@Override
	public Boolean hasEPS()
	{
		return _hasEps;
	}
	
	/**
	 * \brief Compute the active fraction of the bacterium, ignoring EPS
	 * (i.e. only compare active and inert compartments).
	 * 
	 * @return Double noting the fraction of the bacterium that is active.
	 */
	@Override
	public Double getActiveFrac()
	{
		if ( ! _hasInert )
			return 1.0;
		/*
		 * Active fraction is biomass / ( biomass + inert ).
		 */
		int index = this._agentGrid.mySim.getParticleIndex("biomass");
		Double val = particleMass[index];
		index = this._agentGrid.mySim.getParticleIndex("inert");
		val = val / (val + particleMass[index]);
		/*
		 * Check nothing's gone wrong (e.g. if biomass + inert = 0).
		 */
		if ( Double.isNaN(val) )
			val = 1.0;
  		return val;
  	}
	
	/**
	 * \brief Find the particle index for "capsule", i.e. EPS.
	 * 
	 * @return The particle index for "capsule".
	 */
	private int getIndexCapsule()
	{
		return this._agentGrid.mySim.getParticleIndex("capsule");
	}
	
	/**
	 * \brief Send the colour associated to the species to the defined EPS
	 * capsule (if appropriate).
	 * 
	 * @return Color object that this species of Bacterium has been assigned.
	 */
	@Override
	public Color getColorCapsule()
	{
		if ( _epsSpecies == null )
			return getSpeciesParam().epsColor;
		else
			return _epsSpecies.color;
	}
}
