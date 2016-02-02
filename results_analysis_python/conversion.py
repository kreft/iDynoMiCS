#!/usr/bin/python
import sys
import convert_results_for_fig4 as convert

convert.build_life_history(sys.argv[1], 'specific growth rate')
