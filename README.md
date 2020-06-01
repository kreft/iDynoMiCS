# Notes on replicating the results of Wright, Clegg, Coker & Kreft (2020) “Damage repair versus aging in biofilms using an individual-based model”
### [View preprint](https://doi.org/10.1101/2020.01.08.899740)

**We used [iDynoMiCS](https://www.birmingham.ac.uk/generic/idynomics/index.aspx) version 1.5, written in Java, for all simulations, and all analysis scripts were written in Python.**

**The full data for all simulations shown in this manuscript are over 2 TB. We have therefore not provided these data, but have provided protocol files for each type of simulation run in this study. In this repository you can find:**
1. The "biofilm_manuscript_example_protocol" folder contains commented protocol files that explain how these can be modified to produce all protocol files used.
2. The "biofilm_manuscript_results_analysis" folder contains the Python analysis scripts used.
3. All of the iDynoMiCS source code and files necessary to run this (see the above website and associated iDynoMiCS publication - [Lardon et al. (2011)](https://doi.org/10.1111/j.1462-2920.2011.02414.x) - for directions on this).

**Please contact [Robyn Wright](mailto:robyn.wright@dal.ca) with any questions.**

---

Below, we describe the files that we have used here, the main sections that these fall into and how the example protocol files may be edited to reproduce these. All toolbox_ files in the results analysis folder are modified from our previous study, Clegg et al. (2014).

1. **Comparison with Clegg et al. (2014) (where damage accumulation rate, 0.22, is proportional to specific growth rate, which has a maximum of 0.6 h-1)**
	- **Aims of simulations:**
		- Compare the new adaptive repair mechanism with the fixed, optimal repair used by Clegg et al. (2014)
		- Examine the effect of the new adaptive repair mechanism on the growth rates and investment into repair for all six possible repair/division strategies (listed in Figure 1) when growing on their own in the constant environment with no competition (i.e. only one cell is left in the simulation domain).
		- Look at the growth rate distributions for cells of a single strategy at day 100 of simulated growth.
	- **Manuscript:**
		- *Figure 2:* Characteristics of strategies in a constant environment (no competition) showing: (A) investment into repair for the new 
		  adaptive strategy following the old pole cell over many divisions; (B) Specific growth rate of a single cell over consecutive cell divisions 
		  over a 24 h period; and (C) distribution of specific growth rates in populations at steady state (100 days). 
		- *Figure S1:* The investment into repair versus the proportion of protein that has been invested into repair for symmetrically and asymmetrically dividing cells of the adaptive repair strategy.
		- *Figure S2:* Computational replicates of Figure 2. 
		- *Figure S3:* Age and size distributions (scatter plot) for computational replicates of all six strategies at steady state. 
	- **Example protocol files:**
		- follow_old_pole_cell.xml (Figure 2A, 2B, S1, S2)
		- constant_100_day_single.xml (Figure 2C, S1, S2, S3)

	- **Analysis files (and files created):**
		- old_pole_figs_2_S1_S2_S3.py - Files contained in the biofilm_manuscript_results_analysis/following_old_cell folder

2. **Constant and chemostat competitions (where damage accumulation rate, 0.22 h-1, is proportional to specific growth rate, maximum 0.6 h-1; Fig. 3 of the manuscript)**
	- **Aims of simulations:**
		- Determine the fittest strategy for dealing with damage accumulation in spatially homogeneous environments (constant and chemostat)
	- **Manuscript:**
		- *Figure 3:* Time courses showing log biomass ratios of competitions between aging strategies in constant and chemostat environments.
	- **Example protocol files:**
		- constant_competition.xml (Figure 3A-C)
		- chemostat_competition.xml (Figure 3D-F)
	- **Analysis files (and files created):**
		- comps_fig_3.py - Files contained in the biofilm_manuscript_results_analysis/const_chemo_comps folder

3. **Biofilm simulations with no damage accumulation or repair (where specific growth rate is 0.6 h-1 and parameters influencing diffusion and bulk transport are varied)**
	- **Aims of simulations:**
		- Determine parameters that give rise to typical biofilm structures. 
	- **Manuscript:**
		- *Figure S4:* Plot of roughness over different biofilm heights for different values of delta squared. 
		- Plots of biofilms over time are shown in the supplementary [Figshare file](https://figshare.com/articles/Damage_repair_versus_aging_in_biofilms-File_S1_pdf/11520534/1)
	- **Example protocol files:**
		- roughness_test.xml (Figure S4)
		- agent_State_spacing.xml (input file giving initial placement of cells - these can be created using the script even_spacing_to_start.py)
	- **Analysis files (and files created):**
		- roughness_fig_s5_s6.py, roughness_graph.py - Files contained in the biofilm_manuscript_results_analysis/roughness folder

4. **Biofilm competitions (where damage accumulation rate, 0.1 h-1, is not proportional to specific growth rate, maximum 0.6 h-1; Fig. 4 of the manuscript)**
	- **Aims of simulations:**
		- Determine the fittest strategy for cells growing in a biofilm. These were carried out prior to realising that damage accumulation rate must be 
		  made proportional to specific growth rate. 
	- **Manuscript:**
		- *Figure 4:* Damage segregation versus adaptive repair strategies in biofilms with and without ‘styrofoam’.
		- *Figure S5:* Time courses of log biomass ratios for competitions between adaptive repair, fixed repair and damage segregation as well as 
		  control competitions carried out between two cells of the same strategy in biofilms with and without ’styrofoam’. 
		- Plots of all replicate biofilms at the end of the simulations are shown in the supplementary [Figshare file](https://figshare.com/articles/Damage_repair_versus_aging_in_biofilms-File_S1_pdf/11520534/1)
	- **Example protocol file:**
		- styrofoam_competition.xml
	- **Analysis files (and files created):**
		- shrinking_styrofoam_plots_figs_4_supp_file.py - Files contained in the biofilm_manuscript_results_analysis/shrinking_styro_biofilms folder
		- all_stats_analysis_shrinking_styro.py - Files contained in the biofilm_manuscript_results_analysis/stats folder

5. **Determination of which growth rate-proportional damage accumulation rate is equivalent to the previous constant damage accumulation rate (where maximum specific growth rate is set to 0.6 h-1 and damage accumulation rate is either 0.1 h-1 and is not proportional to specific growth rate or is varied between 0-0.25 h-1 and is proportional to specific growth rate)**
	- **Aims of simulations:**
		- Determine which proportional damage accumulation rate is equivalent to a damage accumulation rate of 0.1 h-1 when damage accumulation is not proportional to specific growth rate.
	- **Manuscript:**
		- *Figure S6:* Plots of time in days versus population age, size or growth rate for varying damage accumulation rates.

	- **Analysis files (and files created):**
		- aging_rate_fig_S6.py - Files contained in the biofilm_manuscript_results_analysis/aging_rate folder

6. **Biofilm competitions (where damage accumulation rate, 0.22, is proportional to specific growth rate, maximum 0.6 h-1; Figs. 5-6 of the manuscript)**
	- **Aims of simulations:**
		- Determine the fittest strategy for cells growing in a biofilm where damage accumulation rate is proportional to specific growth rate.
	- **Manuscript:**
		- *Figure 5:* Damage segregation versus adaptive repair strategies in biofilms, showing cells coloured by age and growth rate for one 
		  representative biofilm. 
		- *Figure 6:* Time courses of log biomass ratios between the two competed strategies.
		- *Figure S7:* Time courses of log biomass ratios as well as plots of biofilms for AR/DS (as shown in Figures 5 and 6) and control competitions 
		  of AR/AR and DS/DS.
		- *Table S1:* Results and statistics for all proportional biofilm competition outcomes.
		- Plots of all replicate biofilms at the end of the simulations are shown in the supplementary [Figshare file](https://figshare.com/articles/Damage_repair_versus_aging_in_biofilms-File_S1_pdf/11520534/1)
	- **Example protocol files:**
		- proportional_competition.xml
	- **Analysis files (and files created):**
		- proportional_plots.py - Files contained in the biofilm_manuscript_results_analysis/proportional_biofilms folder
		- combined_biomass_plot_proportional_fig_6.py - biofilm_manuscript_results_analysis/proportional_biofilms/Biomass_time_courses.png
