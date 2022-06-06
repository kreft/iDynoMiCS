#!/usr/bin/python
import sys
import convert_results_for_age_vs_Prep_Ptot as convert


#For Fig4B
#convert.build_population_structure(sys.argv[1], 'specific growth rate')

#For Fig4AC
convert.build_life_history(sys.argv[1], 'age', 'activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair', 'generation')