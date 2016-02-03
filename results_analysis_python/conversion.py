#!/usr/bin/python
import sys
import convert_results_for_fig4 as convert


convert.build_population_structure(sys.argv[1], 'specific growth rate')
