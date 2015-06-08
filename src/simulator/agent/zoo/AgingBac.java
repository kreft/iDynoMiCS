package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.geometry.ContinuousVector;
import utils.ExtraMath;
import utils.LogFile;

/**
 * Laurent Lardon wrote the prototype AgingBac
 * Modified by Jan Kreft (j.kreft@bham.ac.uk) 2009-06-18
 * Modified by Edd Miles
 * Modified by Rob Clegg (rjc096@bham.ac.uk) Oct 2010 - Jan 2011
 * AgingBac is a bacterium with asymmetric or symmetric division driven by parameter alpha
 */
public class AgingBac extends Bacterium
{
	/**
	 * Age is proportion of biomass that is inactive, and has nothing to do
	 * with time!
	 */
	protected Double age = 0.0;
	
	/**
	 * This is different from isDead! We need our own way of dying here,
	 * keeping the dead cell as a ghost that is still present as a particle
	 * and has one active reaction implementing the lysis of the cell
	 * releasing nutrients e.g. COD
	 */
	protected boolean hasDied = false;
	
	/**
	 * 
	 */
	ContinuousVector birthPlace = new ContinuousVector();
	
	/**
	 * 
	 */
	ContinuousVector deathPlace = new ContinuousVector();
	
	/**
	 * 
	 */
	protected Double _deathday = 0.0;
	
	public AgingBac()
	{
		super();
		_speciesParam = new AgingParam();
	}

	public void initFromResultFile(Simulator aSim, String[] singleAgentData)
	{
		// Rewritten by Rob 10/1/11 to agree with subroutines in the agent hierarchy

		// Find the position to start at by using length and number of values read
		// AgingBacv has 9 extra columns: birth (X,Y,Z), age, hasDied, death (X,Y,Z), deathday
		int nValsRead = 9;
		int iDataStart = singleAgentData.length - nValsRead;

		birthPlace.x=Double.parseDouble(singleAgentData[iDataStart]);
		birthPlace.y=Double.parseDouble(singleAgentData[iDataStart+1]);
		birthPlace.z=Double.parseDouble(singleAgentData[iDataStart+2]);
		age         =Double.parseDouble(singleAgentData[iDataStart+3]);
		hasDied     =(Boolean.valueOf(singleAgentData[iDataStart+4])).booleanValue();
		deathPlace.x=Double.parseDouble(singleAgentData[iDataStart+5]);
		deathPlace.y=Double.parseDouble(singleAgentData[iDataStart+6]);
		deathPlace.z=Double.parseDouble(singleAgentData[iDataStart+7]);
		_deathday   =Double.parseDouble(singleAgentData[iDataStart+8]);

		// now go up the hierarchy with the rest of the data
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];

		super.initFromResultFile(aSim,remainingSingleAgentData);
	}

	public void makeKid() throws CloneNotSupportedException
	{
		// Create the new instance
		AgingBac baby = (AgingBac) sendNewAgent();
		
		this._myDivRadius = ExtraMath.deviateFromCV(
				getSpeciesParam().divRadius, getSpeciesParam().divRadiusCV);
		baby._myDivRadius = ExtraMath.deviateFromCV(
				getSpeciesParam().divRadius, getSpeciesParam().divRadiusCV);
		baby._myDeathRadius = ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
				getSpeciesParam().deathRadiusCV);
		
		// Update the lineage
		recordGenealogy(baby);

		// Share mass of all compounds between two daughter cells and compute
		// new size
		divideCompounds(baby, getBabyMassFrac());
		//sonia:chemostat
		if (Simulator.isChemostat){
			// upon division the daughter cells remain with the coordinates of their progenitor

		}else{
			// Compute movement to apply to both cells
			setDivisionDirection(getInteractDistance(baby)/2);

			// move both daughter cells
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}
		// Now register the agent inside the guilds and the agent grid
		baby.registerBirth();
		baby._netVolumeRate = 0.0;
		baby.birthPlace = baby._location;

	}

	public void divideCompounds(AgingBac baby, Double babyMassFrac)
	{
		Double alpha = getSpeciesParam().alpha;
		Double alphaDev = getSpeciesParam().alphaDev;
		/*
		 * If there's any stochasticity in segregation, determine alpha now.
		 */
		if ( alphaDev > 0.0 )
		{
		    Double alphaTemp;
		    do {
			alphaTemp = alpha + alphaDev*ExtraMath.getNormRand();
		    } while ( Math.abs(alphaTemp) > 1.0 );
		    alpha = alphaTemp;
		}
		/*
		 * Need to store pre-division mother's protein levels so they can be
		 * correctly distributed between daughter cells. This will then ensure
		 * a mass balance in undamaged, damaged and total protein masses.
		 */
		Double totalDam = this.particleMass[1];
		Double totalUndam = this.particleMass[0];
		Double totalProt = totalDam + totalUndam;
		/*
		 *  babyMassFrac = R, mumMassFrac = (1-R)
		 *  damFactor = DRa
		 */
		Double mumMassFrac = 1 - babyMassFrac;
		Double damFactor = totalDam * babyMassFrac * alpha;
		/*
		 * Damage segregation equations (modified) from Erjavec et al (2008)
		 * (equations 5 & 6). Note that 0<alpha<1 where 0 means symmetry,
		 * and 1 means complete asymmetry.
		 * 
		 * Modification is dealing with a post-division mother full of damage.
		 * Warning note to log that an overflow has occurred. You can search
		 * for this in the log0.txt file using the keyword "overflow" 
		 */
		if ( totalDam * alpha * (1-babyMassFrac) > totalUndam * babyMassFrac )
			LogFile.writeLog("Warning: overflow of damage into younger sibling has occurred.");
		/*
		 * Baby's damaged = max[ DR-DRa, D-(1-R)T -> any overflow damage]
		 */
		baby.particleMass[1] = Math.max(totalDam*babyMassFrac-damFactor,
											totalDam-mumMassFrac*totalProt); 
		/*
		 * Baby's undamaged = min[ UR +DRa, U -> all undamaged]
		 */
		baby.particleMass[0] = Math.min(babyMassFrac * totalUndam + damFactor,
																totalUndam);
		/*
		 * Parent's damaged = min[ D(1-R)+DRa, RT -> all possible damage]
		 */
		this.particleMass[1] = Math.min(totalDam * mumMassFrac + damFactor,
													mumMassFrac * totalProt);
		/*
		 * Parent's undamaged = max[ U(1-R)-DRa, 0 -> no space for undamaged]
		 */
		this.particleMass[0] = Math.max(totalUndam * mumMassFrac - damFactor,
																		0.0);
		/*
		 * If alpha has deviated to a negative value, we need to check
		 * everything's still ok.
		 */
		//  
		if ((totalUndam + alpha*totalDam < 0)||((1 - alpha)*babyMassFrac > 1))
		{
		    baby.particleMass[1] = Math.min(totalDam, babyMassFrac*totalProt);
		    baby.particleMass[0] = Math.max(0.0,
		    							babyMassFrac * totalProt - totalDam);
		    this.particleMass[1] = Math.max(0.0,
		    							totalDam - babyMassFrac * totalProt);
		    this.particleMass[0] = Math.min(mumMassFrac*totalProt, totalUndam);
		}
		/*
		 * Modified from Edd's final files
		 */
		if ( this.particleMass[0] < 0.0 || baby.particleMass[0] < 0.0 ||
				this.particleMass[1] < 0.0 || baby.particleMass[1] < 0.0 )
		{
		    System.out.println(" ");
		    System.out.println("Error:");
		    System.out.println("  predivision cell's undamaged "+totalUndam+
		    											", damaged "+totalDam);
		    System.out.println("  parent's undamaged "+this.particleMass[0]+
		    								", damaged "+this.particleMass[1]);
		    System.out.println("  baby's undamaged "+baby.particleMass[0]+
		    								", damaged "+baby.particleMass[1]);
		    System.out.println("alpha "+alpha+" babyMassFrac "+babyMassFrac);
		    System.out.println(" ");
		    System.exit(-1);
		}
		if (Math.abs(this.particleMass[1]+baby.particleMass[1]-totalDam) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's damaged ("+this.particleMass[1]+")");
			System.out.println(" + baby's damaged ("+baby.particleMass[1]+") != mother's damaged ("+totalDam+")");
			System.exit(-1);
		}
		if(Math.abs(this.particleMass[0]+baby.particleMass[0]-totalUndam) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's undamaged ("+this.particleMass[1]+")");
			System.out.println(" + baby's undamaged ("+baby.particleMass[1]+") != mother's undamaged ("+totalDam+")");
			System.exit(-1);
		}	
		/*
		 * Age 0 is youngest, 1 oldest, fraction of inactive biomass
		 * (particleMass[1]) over total biomass.
		 */
		this.updateAge(); 
		baby.updateAge();
		/*
		 * Update radii, masses and volumes.
		 */
		this.updateSize();
		baby.updateSize();
		this.updateGrowthRates();
		baby.updateGrowthRates();
	}
	
	public void updateAge()
	{
		if ( this._totalMass < 0.0 )
		{
			LogFile.writeLogAlways("Warning: _totalMass "+this._totalMass+
										", my id"+_family+", "+_genealogy);
			die(true);
			return;
		}
		if ( this.particleMass[0] < 0.0 )
		{
			LogFile.writeLogAlways("Warning: particleMass[0]"+
					this.particleMass[0]+", my id"+_family+", "+_genealogy);
			die(true);
			return;
		}
		age = particleMass[1] / (particleMass[0]+particleMass[1]);
	}


	/**
	 * This is a really important override of willDie().
	 * Apart from the instances where an AgingBac's mass is negative (i.e. something's 
	 * gone wrong!), it tells any AgingBac whose age (not chronological, but proportion 
	 * of total mass that is inactive) is more than 99.9% to switch on hasDied, convert 
	 * all remaining active biomass to inactive, set its deathPlace and deathday, 
	 * unregister from all active reactions and begin lysing.
	 */
	@Override
	public boolean willDie()
	{
		/*
		 * We only die because of old age, no other causes of death should
		 * muddy the waters.
		 */
		if ( this._totalMass < 0.0 )
		{
			LogFile.writeLogAlways("Warning: _totalMass "+this._totalMass+
											", my id"+_family+", "+_genealogy);
			return true;
		}
		if ( this.particleMass[0] < 0.0 )
		{
			LogFile.writeLogAlways("Warning: particleMass[0]"+
					this.particleMass[0]+", my id"+_family+", "+_genealogy);
			return true;
		}
		/*
		 * If hasDied, we have done everything already.
		 */
		if ( hasDied )
			return false; 
		/*
		 * Check if at least 99.9% damaged - if so, round up to 100% 
		 */
		if (age >= 0.999)
		{
			hasDied = true;
			/*
			 * We convert all ~0.001 remaining active mass to inactive.
			 */
			this.particleMass[1] += this.particleMass[0];
			this.particleMass[0] = 0.0;
			this.age = 1.0;
			/*
			 * Make note of _deathDay and deathPlace.
			 */
			this.deathPlace = this._location;
			this._deathday = SimTimer.getCurrentTime();
			LogFile.writeLog("DIED! age "+age+", id "+_family+", "+_genealogy+
															", "+_generation);
			/*
			 * Turn off all active reactions.
			 */
			unregisterFromAllActiveReactions();
			/*
			 * Turn on the lysis of dead cells reaction, if this reaction exists.
			 */
			if (_agentGrid.mySim.getReactionIndex("lysis") > -1)
				switchOnReaction(_agentGrid.mySim.getReaction("lysis"));
		}
		return false;
	}
	
	@Override
	public void grow()
	{
		updateGrowthRates();
		for (int i = 0; i < particleMass.length; i++)
			particleMass[i] += deltaParticle[i];
	}
	
	@Override
	protected void updateGrowthRates()
	{
		super.updateGrowthRates();
		repair();
	}
	
	protected void repair()
	{
		Double beta = getSpeciesParam().beta;
		Double pDam = this.particleMass[1];
		if ( beta == 0.0 || pDam == 0.0 )
			return;
		/*
		 * Apply the changes in active and damaged biomass. Note that here it
		 * is assumed that active and damaged biomass have the same densities.
		 */
		Double tStep = SimTimer.getCurrentTimeStep();
		Double rMax = getSpeciesParam().rMax;
		Double pAct =  this.particleMass[0];
		Double rate = rMax * pDam / (beta * pAct + pDam);
		Double repMass = beta * pAct * Math.expm1(tStep * rate);
		Double repY = getSpeciesParam().repY;
		deltaParticle[0] += repY * repMass;
		deltaParticle[1] -= repMass;
		_netGrowthRate -= (1 - repY) * beta * pAct * rate;
		_netVolumeRate -= (1 - repY) * beta * pAct * rate /
										getSpeciesParam().particleDensity[0];
	}
	
	/**
	 * Called at each time step (under the control of the method Step of the
	 * class Agent to avoid multiple calls)
	 */
	@Override
	protected void internalStep()
	{
		/*
		 * Compute mass growth over all compartments
		 */
		grow();
		updateSize();
		/*
		 * This is specific to AgingBac
		 */
		if ( ! hasDied )
			updateAge();
		/*
		 * Test if the EPS capsule has to be excreted
		 */
		manageEPS();
		/*
		 * Divide if you have to EDD:
		 */
		if ( willDivide() && ! hasDied )
			divide();
		/*
		 * Die if you have to. Here we kill off only those whose mass is
		 * non-positive (see first 2 cases in willDie() and the catch in
		 * updateAge()).
		 */
		if ( willDie() )
			die(true);
	}
	
	/**
	 * This will be added to the END of the string from higher classes.
	 * Edd added death & Rob (17/1/11) added deathDay.
	 */
	public StringBuffer sendHeader()
	{
		StringBuffer tempString = new StringBuffer(super.sendHeader());
		tempString.append(",birthX,birthY,birthZ,age,hasDied");
		tempString.append(",deathX,deathY,deathZ,deathday");
		return tempString;
	}
	
	/**
	 * This will be added to the END of the string from higher classes.
	 * Edd added death & Rob (17/1/11) added deathDay.
	 */
	public StringBuffer writeOutput()
	{
		StringBuffer tempString = new StringBuffer(super.writeOutput());
		tempString.append(","+birthPlace.x+","+birthPlace.y+","+birthPlace.z+
					","+age+","+hasDied+","+deathPlace.x+","+deathPlace.y+","+
					deathPlace.z+","+_deathday);
		return tempString;
	}	

	public AgingParam getSpeciesParam()
	{
		return (AgingParam) _speciesParam;
	}
}
