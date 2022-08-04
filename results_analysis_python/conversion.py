#!/usr/bin/python
import sys
import convert_results_for_fig4 as convert


#For Fig4B
#convert.build_population_structure(sys.argv[1], 'specific growth rate')

#For Fig4AC
convert.build_life_history(sys.argv[1], 'specific growth rate', 'age')
