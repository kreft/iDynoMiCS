#!/usr/bin/python
from __future__ import division
#import aging_extras
import os
#import math
#from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results
#from mpl_toolkits.axes_grid1 import make_axes_locatable

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'comparison_clegg', 'Final_versions', 'Reproducing_figure_4', 'nontoxic')
output_path = os.path.join(input_path, 'figure_4_life_history_B_NT.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')
fs = 8 

axisB = fig.add_subplot('B', 111)
paths = ['NT010AS12NR100_results.xml', 'NT010AS12007100_results.xml', 'NT010AS12OB100_results.xml']
attributes = {'name':'age population structure',
                    'starting_time':'2400', 'header':'bin,frequency'}
colorsB = ['red','#F781F3','#fbbf07']

for i in range(3):
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        bin_edges = result.vars['bin'].split('-')
        result.vars['mid_point'] = (float(bin_edges[0])+float(bin_edges[1]))/2
    result_set.members.sort(key=lambda r: r.vars['mid_point'])
    x = [float(r.vars['frequency']) for r in result_set.members]
    y = [r.vars['mid_point'] for r in result_set.members]
    axisB.plot(x, y)

'''
axisB.text(0.030,0.57,'0', color=colorsB[2], fontsize=fs)
axisB.text(0.044,0.51,'1', color=colorsB[2], fontsize=fs)
axisB.text(0.080,0.38,'2', color=colorsB[2], fontsize=fs)
axisB.text(0.025,0.27,'3-', color=colorsB[2], fontsize=fs)
axisB.text(0.026,0.57,'0', color=colorsB[0], fontsize=fs)
axisB.text(0.040,0.51,'1', color=colorsB[0], fontsize=fs)
axisB.text(0.075,0.38,'2', color=colorsB[0], fontsize=fs)
axisB.text(0.028,0.24,'3', color=colorsB[0], fontsize=fs)
axisB.text(0.014,0.08,'4', color=colorsB[0], fontsize=fs)
axisB.text(0.016,0.00,'5', color=colorsB[0], fontsize=fs)
axisB.text(0.045,0.54,'0', color=colorsB[1], fontsize=fs)
axisB.text(0.055,0.48,'1', color=colorsB[1], fontsize=fs)
axisB.text(0.095,0.40,'2', color=colorsB[1], fontsize=fs)
axisB.text(0.024,0.33,'3', color=colorsB[1], fontsize=fs)


aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
'''
#axisB.set_xlim([0, 0.105])
#axisB.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10])
#axisB.set_ylim([0, 0.6])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate Distribution')
#Only needed when no A axis - also make y tick labels true if A absent
#axisB.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
#setp( axisB.get_yticklabels(), visible=False)

fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)