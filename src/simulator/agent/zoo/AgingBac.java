package simulator.agent.zoo;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SpatialGrid;
import simulator.geometry.ContinuousVector;
import utils.ExtraMath;
import utils.LogFile;

/**
 * There are a couple of different ways to call the reactions , depending on what you want.
 * reactionActive is a list of all of the active reactions, compiled from this. So you can
 * either index the actual reaction that you want, or index the active reaction you want.
 * Indexing starts at 0, and is [r][c]
 * First two reactions in protocols are two different growth strategies and are called:
 * 				allReactions[reactionActive.get(0)]
 * 			growthN allReactions[0] 
 * 			growthR allReactions[1]
 * 			Growth can only be one or the other of these, so both are active reaction 0.
 * Third reaction is aging rate, so could be called:
 * 				allReactions[2] or allReactions[reactionActive.get(1)]
 * Fourth reaction is removal of active biomass - this was for comparison with Erjavec and
 * default value is false:
 * 				allReactions[3] or allReactions[reactionActive.get(2)]
 * Fifth reaction is removal of damaged biomass - again Erjavec and default false:
 * 				allReactions[4] or allReactions[reactionActive.get(3)]
 * 
 * I will stick with using the active reaction index where possible. 
 * 
 * Now to include 4 different types of biomass: active biomass invested in growth, active
 * biomass invested in repair, inactive biomass invested in growth and inactive biomass
 * invested in repair.
 * 
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
		 * a mass balance in undamaged growth, undamaged repair, damaged growth, damaged repair and total protein masses.
		 */
		Double totalUndamGrowth = this.particleMass[0];
		Double totalUndamRepair = this.particleMass[1];
		Double totalDamGrowth = this.particleMass[2];
		Double totalDamRepair = this.particleMass[3];
		Double totalDam = totalDamGrowth+totalDamRepair;
		Double totalUndam = totalUndamGrowth+totalUndamRepair;
		Double totalGrowth = totalUndamGrowth+totalDamGrowth;
		Double totalRepair = totalUndamRepair+totalDamRepair;
		Double totalProt = totalDam+totalUndam;
		/*
		 *  babyMassFrac = R, mumMassFrac = (1-R)
		 *  damFactor = DRa
		 */
		Double mumMassFrac = 1 - babyMassFrac;
		Double damFactorGrowth = totalDamGrowth * babyMassFrac * alpha;
		Double damFactorRepair = totalDamRepair * babyMassFrac * alpha;
		/*
		 * Damage segregation equations (modified) from Erjavec et al (2008)
		 * (equations 5 & 6). Note that 0<alpha<1 where 0 means symmetry,
		 * and 1 means complete asymmetry.
		 *
		 * Modification is dealing with a post-division mother full of damage
		 * This is to make sure that we have no negative values. If running with stochasticity
		 * then you may need to check that none of the other terms are negative also.
		 */
		if ( mumMassFrac*totalUndamGrowth < damFactorGrowth || 
								mumMassFrac*totalUndamRepair < damFactorRepair)
		{
			/* Warning note to log that an overflow has occurred. You can search
			 * for this in the log0.txt file using the keyword "overflow" */
			LogFile.writeLog("Warning: overflow of damage into younger sibling has occurred.");
			
			/* Baby's damaged growth = D-(1-R)T -> any overflow damage */
			baby.particleMass[2] = totalDamGrowth - mumMassFrac*totalGrowth;
			/* Baby's damaged repair = D-(1-R)T -> any overflow damage */
			baby.particleMass[3] = totalDamRepair - mumMassFrac*totalRepair;
			/* Baby's undamaged growth = U -> all undamaged */
			baby.particleMass[0] = totalUndamGrowth;
			/* Baby's undamaged repair = U -> all undamaged */
			baby.particleMass[1] = totalUndamRepair;
			/* Parents damaged growth = RT -> all possible damage */
			this.particleMass[2] = mumMassFrac*totalGrowth;
			/* Parents damaged repair = RT -> all possible damage */
			this.particleMass[3] = mumMassFrac*totalRepair;
			/* Parents undamaged growth = 0 -> no space for undamaged */
			this.particleMass[0] = 0.0;
			/* Parents undamaged repair = 0 -> no space for undamaged */
			this.particleMass[1] = 0.0;
		}
		else
		{
			/* Baby's damaged growth = DR-DRa */
			baby.particleMass[2] = totalDamGrowth*babyMassFrac - damFactorGrowth;
			/* Baby's damaged repair = DR-DRa */
			baby.particleMass[3] = totalDamRepair*babyMassFrac - damFactorRepair;
			/* Baby's undamaged growth = UR +DRa */
			baby.particleMass[0] = babyMassFrac*totalUndamGrowth + damFactorGrowth;
			/* Baby's undamaged repair = UR +DRa */
			baby.particleMass[1] = babyMassFrac*totalUndamRepair + damFactorRepair;
			/* Parents damaged growth = D(1-R)+DRa */
			this.particleMass[2] = totalDamGrowth*mumMassFrac + damFactorGrowth;
			/* Parents damaged repair = D(1-R)+DRa */
			this.particleMass[3] = totalDamRepair*mumMassFrac + damFactorRepair;
			/* Parents undamaged growth = U(1-R)-DRa */
			this.particleMass[0] = totalUndamGrowth*mumMassFrac - damFactorGrowth;
			/* Parents undamaged repair = U(1-R)-DRa */
			this.particleMass[1] = totalUndamRepair*mumMassFrac - damFactorRepair;
		}
		
		/*
		 * If alpha has deviated to a negative value, we need to check
		 * everything's still ok.
		 */
		//  
		if ((totalUndam + alpha*totalDam < 0) || ((1 - alpha)*babyMassFrac > 1))
		{
			baby.particleMass[2] = Math.min(totalDamGrowth, babyMassFrac*totalGrowth);
			baby.particleMass[3] = Math.min(totalDamRepair, babyMassFrac*totalRepair);
			baby.particleMass[0] = Math.max(0.0, babyMassFrac*totalGrowth-totalDamGrowth);
			baby.particleMass[1] = Math.max(0.0, babyMassFrac*totalRepair-totalDamRepair);
			this.particleMass[2] = Math.max(0.0, totalDamGrowth-babyMassFrac*totalGrowth);
			this.particleMass[3] = Math.max(0.0, totalDamRepair-babyMassFrac*totalRepair);
			this.particleMass[0] = Math.min(mumMassFrac*totalGrowth, totalUndamGrowth);
			this.particleMass[1] = Math.min(mumMassFrac*totalRepair, totalUndamRepair);
		}
		/*
		 * Modified from Edd's final files
		 */
		if ( this.particleMass[0] < 0.0 || baby.particleMass[0] < 0.0 ||
				this.particleMass[1] < 0.0 || baby.particleMass[1] < 0.0 ||
				this.particleMass[2] < 0.0 || baby.particleMass[2] < 0.0 ||
				this.particleMass[3] < 0.0 || baby.particleMass[3] < 0.0)
		{
		    System.out.println(" ");
		    System.out.println("Error:");
		    System.out.println("  predivision cell's undamaged growth "+totalUndamGrowth+
		    											", damaged growth "+totalDamGrowth);
		    System.out.println("  predivision cell's undamaged repair "+totalUndamRepair+
					", damaged repair "+totalDamRepair);
		    
		    System.out.println("  parent's undamaged growth "+this.particleMass[0]+
		    								", damaged growth "+this.particleMass[2]);
		    System.out.println("  parent's undamaged repair "+this.particleMass[1]+
					", damaged repair "+this.particleMass[3]);
		    
		    System.out.println("  baby's undamaged growth "+baby.particleMass[0]+
		    								", damaged growth "+baby.particleMass[2]);
		    System.out.println("  baby's undamaged repair "+baby.particleMass[1]+
					", damaged repair "+baby.particleMass[3]);
		    
		    System.out.println("alpha "+alpha+" babyMassFrac "+babyMassFrac);
		    System.out.println(" ");
		    System.exit(-1);
		}
		if (Math.abs(this.particleMass[2]+baby.particleMass[2]-totalDamGrowth) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's damaged growth ("+this.particleMass[2]+")");
			System.out.println(" + baby's damaged growth ("+baby.particleMass[2]+") != mother's damaged growth ("+totalDamGrowth+")");
			System.exit(-1);
		}
		if (Math.abs(this.particleMass[3]+baby.particleMass[3]-totalDamRepair) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's damaged repair ("+this.particleMass[3]+")");
			System.out.println(" + baby's damaged repair ("+baby.particleMass[3]+") != mother's damaged repair ("+totalDamRepair+")");
			System.exit(-1);
		}
		if(Math.abs(this.particleMass[0]+baby.particleMass[0]-totalUndamGrowth) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's undamaged growth ("+this.particleMass[0]+")");
			System.out.println(" + baby's undamaged growth ("+baby.particleMass[0]+") != mother's undamaged growth ("+totalUndamGrowth+")");
			System.exit(-1);
		}
		if(Math.abs(this.particleMass[1]+baby.particleMass[1]-totalUndamRepair) > 1E-9)
		{
			System.out.println("Error: Incorrect protein levels post partitioning! Parent's undamaged repair ("+this.particleMass[1]+")");
			System.out.println(" + baby's undamaged repair ("+baby.particleMass[1]+") != mother's undamaged repair ("+totalUndamRepair+")");
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
		if ( this.particleMass[0]+this.particleMass[1] < 0.0 )
		{
			LogFile.writeLogAlways("Warning: particleMass[0]"+
					this.particleMass[0]+"+ particleMass[1]"+this.particleMass[1]+", my id"+_family+", "+_genealogy);
			die(true);
			return;
		}
		age = (particleMass[2]+particleMass[3]) / 
			 (particleMass[0]+particleMass[1]+particleMass[2]+particleMass[3]);
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
		if ( this.particleMass[1] < 0.0 )
		{
			LogFile.writeLogAlways("Warning: particleMass[1]"+
					this.particleMass[1]+", my id"+_family+", "+_genealogy);
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
			this.particleMass[2] += this.particleMass[0];
			this.particleMass[3] += this.particleMass[1];
			this.particleMass[0] = 0.0;
			this.particleMass[1] = 0.0;
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
		for (int i = 0; i < particleMass.length; i++)
			deltaParticle[i] = 0.0;
		
		double tStep = SimTimer.getCurrentTimeStep();
		double pActGrowth =  this.particleMass[0];
		double pActRepair = this.particleMass[1];
		double pDamGrowth = this.particleMass[2];
		double pDamRepair = this.particleMass[3];
		double age = this.age;
		double Mu = allReactions[reactionActive.get(0)].computeSpecGrowthRate(this);
		double a = allReactions[reactionActive.get(1)].getKinetic()[0];
		double beta;
		if ( this.getSpeciesParam().isOptimalRepair )
			beta = optB();
		else
			beta = this.getSpeciesParam().beta;
		/*
		 * Growth
		 */
		/* First work out the net rates (independent of tStep). */
		double temp = Mu;
		if ( this.getSpeciesParam().isToxic )
			temp *= ( 1 - age );
		_netGrowthRate = temp;
		if ( ! this.getSpeciesParam().isLinear )
			_netGrowthRate *= pActGrowth;
		_netVolumeRate = _netGrowthRate/getSpeciesParam().particleDensity[0];
		/* Now work out the deltaParticles and netVolumeRate. */
		temp *= tStep;
		if ( ! this.getSpeciesParam().isLinear )
			temp = pActGrowth * Math.expm1(temp);
		deltaParticle[0] += (1-beta) * temp;
		deltaParticle[1] += beta * temp;
		/*
		 * Aging & repair
		 */
		Double ageMass = pActGrowth * Math.expm1(-a*tStep);
		deltaParticle[0] += ageMass;
		deltaParticle[2] -= ageMass;
		ageMass = pActRepair * Math.expm1(-a*tStep);
		deltaParticle[1] += ageMass;
		deltaParticle[3] -= ageMass;
		repair();
		/*
		 * Removal (Erjavec comparison)
		 */
		if (allReactions.length > 3)
		{
			Double k2 = allReactions[3].getKinetic()[0];
			Double k4 = allReactions[4].getKinetic()[0];
			deltaParticle[0] += pActGrowth * Math.expm1( - k2 * tStep );
			deltaParticle[1] += pActRepair * Math.expm1( - k2 * tStep );
			deltaParticle[2] += pDamGrowth * Math.expm1( - k4 * tStep );
			deltaParticle[3] += pDamRepair * Math.expm1(- k4  * tStep );
			_netGrowthRate -= pActGrowth * k2;
			_netGrowthRate -= pActRepair * k2;
			_netGrowthRate -= pDamGrowth * k4;
			_netGrowthRate -= pDamRepair * k4;
		}
		_netVolumeRate = _netGrowthRate/(getSpeciesParam().particleDensity[0]);
	}

	
	protected double optB()
	{
		double repY = getSpeciesParam().repY;
		double Mu = allReactions[reactionActive.get(0)].computeSpecGrowthRate(this);
		double optB = repY / Mu;
		if ( this.getSpeciesParam().isToxic )
			optB *= 1/(1-age);
		return (age/(1-age)) * ( Math.sqrt(optB) - 1);
	}
	
	
	protected void repair()
	{
		Double pActRepair = this.particleMass[1];
		Double pDamGrowth = this.particleMass[2];
		Double pDamRepair = this.particleMass[3];
		Double pDam = pDamGrowth + pDamRepair;
		if ( pActRepair == 0.0 || pDam == 0.0 )
			return;
		/*
		 * Apply the changes in active and damaged biomass. Note that here it
		 * is assumed that active and damaged biomass have the same densities.
		 */
		Double tStep = SimTimer.getCurrentTimeStep();
		Double rMax = getSpeciesParam().rMax;
		Double specRate = rMax / ( pActRepair + pDam );
		Double rateGrowth = pDamGrowth * specRate;
		Double rateRepair = pDamRepair * specRate;
		Double repMassGrowth = pActRepair * Math.expm1(tStep * rateGrowth);
		Double repMassRepair = pActRepair * Math.expm1(tStep * rateRepair);
		Double repY = getSpeciesParam().repY;
		deltaParticle[0] += repY * repMassGrowth;
		deltaParticle[1] += repY * repMassRepair;
		deltaParticle[2] -= repMassGrowth;
		deltaParticle[3] -= repMassRepair;
		_netGrowthRate -= (1 - repY) * pActRepair * rateGrowth;
		_netGrowthRate -= (1 - repY) * pActRepair * rateRepair;
		_netVolumeRate -= (1 - repY) * pActRepair * (rateGrowth + rateRepair) /
										getSpeciesParam().particleDensity[0];
	}
	
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex)
	{
		if ( isDead )
			return;
		Double value = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		if ( this.getSpeciesParam().isToxic )
			value *= (1-this.age);
		aSpG.addValueAt(value, _location);
	}
	
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG)
	{
		if ( isDead )
			return;
		Double value = _totalMass/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		if ( this.getSpeciesParam().isToxic )
			value *= (1-this.age);
		aSpG.addValueAt(value, _location);
	}
	
	@Override
	public boolean willDivide()
	{
		// this ensures that the checks for when to divide don't occur too often;
		// at most they will occur at the rate of AGENTTIMESTEP
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if ( _timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep() )
			return false;

		// at this point we will actually check whether to divide
		_timeSinceLastDivisionCheck = 0.0;
		
		double checkRad;
		
		if ( this.getSpeciesParam().isPint )
		{
			checkRad = ExtraMath.radiusOfASphere(particleMass[0] /
										getSpeciesParam().particleDensity[0]);
		}
		else
		{
			checkRad = getRadius(false);
		}
		return checkRad  > this._myDivRadius;
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
