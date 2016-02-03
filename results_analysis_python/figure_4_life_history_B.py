#!/usr/bin/python
from __future__ import division
import aging_extras
import os
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results

#base_path = os.path.join('~', 'Dropbox', 'EclipseWorkspace', 'iDynoAge')
base_path = os.path.join('C:\\', 'Users', 'RJW598', 'git', 'iDynoMiCS')
#base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results')
output_path = os.path.join(base_path, 'figure_4_life_history.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')

axisB = fig.add_subplot('', 111)
paths = ['an10_const_steady_results.xml', 'ao10_const_steady_results.xml', 'T010AS12OBConstEnvH110D_results.xml']
attributes = {'name':'specific growth rate population structure',
                    'starting_time':'2400', 'header':'bin,frequency'}
colors = ['red','#F781F3','#fbbf07']
fs = 8

xlim, ylim = 0.082, 0.6

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
    axisB.plot(x, y, color=colors[i])
axisB.text(0.021, 0.575,'0', color=colors[0], fontsize=fs)
axisB.text(0.055, 0.36, '1', color=colors[0], fontsize=fs)
axisB.text(0.019, 0.245, '2', color=colors[0], fontsize=fs)
axisB.text(0.009, 0.08, '3', color=colors[0], fontsize=fs)
axisB.text(0.009, 0.0, '4', color=colors[0], fontsize=fs)
axisB.text(0.03, 0.515,'0', color=colors[1], fontsize=fs)
axisB.text(0.075, 0.38, '1', color=colors[1], fontsize=fs)
axisB.text(0.012, 0.35, '2', color=colors[1], fontsize=fs)
axisB.text(0.0155, 0.31, '3', color=colors[1], fontsize=fs)
axisB.text(0.006, 0.24, '4', color=colors[1], fontsize=fs)
axisB.text(0.003, 0.16, '5', color=colors[1], fontsize=fs)
aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
axisB.set_xlim([0.00, xlim])
axisB.set_xticks([0.00, 0.02, 0.04, 0.06, 0.08])
axisB.set_ylim([0, ylim])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate Distribution')
axisB.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
setp( axisB.get_yticklabels(), visible=False)

fig.process_subplots()
#axisB.tick_params(left='off')
#axisD.tick_params(left='off')
#axisCD.tick_params(top="off", bottom="off", right="off", left="off")
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)