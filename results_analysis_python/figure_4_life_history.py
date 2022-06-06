#!/usr/bin/python
from __future__ import division
import aging_extras
import os
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results

base_path = os.path.join('~', 'Dropbox', 'EclipseWorkspace', 'iDynoAge')
base_path = toolbox_basic.check_path(base_path)
input_path = os.path.join(base_path, 'results_files', 'figure_4')
output_path = os.path.join(base_path, 'figures', 'figure_4_life_history.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')

axisA = fig.add_subplot('A', 221)
paths = ['an10_life_history_results.xml', 'ao10_life_history_results.xml',
         'sn10_life_history_results.xml', 'so10_life_history_results.xml']
colors = ['red','#F781F3','blue','#81BEF7']
attributes = {'name':'specific growth rate', 'header':'time,value'}

xlim, ylim = 24, 0.6

for i in range(4):
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    val = [float(r.vars['value']) for r in result_set.members]
    j = 0
    while j < len(t)-1:
        x = [t[j]]
        y = [val[j]]
        while (j < len(t)-1) and (abs(val[j+1]-val[j]) < 0.05):
            x.append(t[j+1])
            y.append(val[j+1])
            j += 1
        axisA.plot(x, y, color=colors[i])
        j += 1

fs = 8 

axisA.text(9, 0.36, 'Symmetric, No repair', color=colors[2], va='top',
                                                                   fontsize=fs)
axisA.text(9, 0.42, 'Symmetric, Optimal repair', color=colors[3], va='top',
                                                                   fontsize=fs)
axisA.text(1,0.56,'Asymmetric, No repair', color=colors[0], va='center',
                                                                   fontsize=fs)
axisA.text(2,0.48,'Asymmetric, Optimal repair', color=colors[1], va='center',
                                                                   fontsize=fs)
axisA.text(1.5,0.51,'0', fontsize=fs)
axisA.text(3.2, 0.395, '1', fontsize=fs)
axisA.text(3,0.25,'2',color=colors[0], fontsize=fs)
axisA.text(6,0.1,'3',color=colors[0], fontsize=fs)
axisA.text(4.5,0.32,'2',color=colors[1], fontsize=fs)
axisA.text(7,0.26,'3',color=colors[1], fontsize=fs)
axisA.text(10,0.21,'4',color=colors[1], fontsize=fs)
axisA.text(14,0.16,'5',color=colors[1], fontsize=fs)
axisA.text(18,0.11,'6',color=colors[1], fontsize=fs)
axisA.text(22,0.03,'7',color=colors[1], fontsize=fs)
aging_extras.draw_cell(axisA, 0.59*xlim, 0.69*xlim, 0.9*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisA, 0.79*xlim, 0.8*ylim, 0.99*xlim, ylim, y2xscale=ylim/xlim)
axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 4, 8, 12, 16, 20])#, 24])
axisA.set_ylim([0, ylim])
axisA.set_xlabel(r'Time (h)')
axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
axisA.set_title('Following Old Pole Cell Over Time')


axisB = fig.add_subplot('B', 222)
paths = ['an10_const_steady_results.xml', 'ao10_const_steady_results.xml']
attributes = {'name':'specific growth rate population structure',
                    'starting_time':'2400', 'header':'bin,frequency'}

xlim, ylim = 0.082, 0.6

for i in range(2):
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
setp( axisB.get_yticklabels(), visible=False)


axisCD = fig.add_subplot('', 212, frameon=False)

ms = 2
axisC = fig.add_subplot('C', 223)
paths = ['an10_const_toxic_agent_state_last.xml',
         'sn10_const_toxic_agent_state_last.xml']
colorsC = [colors[0], colors[2]]
for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    agent_output = toolbox_results.AgentOutput(path=results_path)
    species = toolbox_results.SpeciesOutput(agent_output)
    aging_extras.scatter_population(axisC, species, 'totalBiomass', 'age',
                                                color=colorsC[i], markersize=ms)
width, height = 5, 0.12
left, right, bottom = 285, 627, -0.01
'''
axisC.plot([left + width, left, left, left + width],
                                          [bottom]*2 + [bottom+height]*2, 'k-')
axisC.plot([right - width, right, right, right - width],
                                          [bottom]*2 + [bottom+height]*2, 'k-')
'''
axisC.text(290, 0.08, 'Generation 0',
                            color=colors[0], fontsize=fs, ha='left', va='top')
axisC.text(290, 0.18, '1', color=colors[0], fontsize=fs, ha='left', va='top')
axisC.text(290, 0.37, '2', color=colors[0], fontsize=fs, ha='left', va='top')
axisC.text(290, 0.6, '3', color=colors[0], fontsize=fs, ha='left', va='top')
axisC.text(290, 0.95, '4', color=colors[0], fontsize=fs, ha='left', va='top')
aging_extras.draw_cell(axisC, 488, 521, 0.9, 16.7, y2xscale=1/333, toxic=True,
                                                           arrow=3.3, repair=0)
aging_extras.draw_const_env(axisC, 554, 0.8, 621, 1.0, y2xscale=1/333)
axisC.set_xlim([285, 630])
axisC.set_ylim([0.0, 1.0])
axisC.set_title('Age & Size Distributions')
axisC.set_ylabel(r'Cellular age $P_{dam}/(P_{act}+P_{dam})$')
# Legend
left, right = 451, 621
top, bottom = 0.78, 0.62
asym_y = 0.735
sym_y = 0.665
axisC.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
axisC.plot([left + 10], [asym_y], 'o', color=colorsC[0], markeredgecolor='none', markersize=ms)
axisC.text(left + 20, asym_y, 'Asymmetric, No repair', va='center', ha='left', fontsize=fs)
axisC.plot([left + 10], [sym_y], 'o', color=colorsC[1], markeredgecolor='none', markersize=ms)
axisC.text(left + 20, sym_y, 'Symmetric, No repair', va='center', ha='left', fontsize=fs)


axisD = fig.add_subplot('D', 224)
paths = ['ao10_const_toxic_agent_state_last.xml',
         'so10_const_toxic_agent_state_last.xml']
colorsD = [colors[1], colors[3]]
for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    agent_output = toolbox_results.AgentOutput(path=results_path)
    species = toolbox_results.SpeciesOutput(agent_output)
    aging_extras.scatter_population(axisD, species, 'totalBiomass', 'age', color=colorsD[i])
axisD.text(285, 0.05, '0', color=colors[1], fontsize=fs, ha='left', va='top')
axisD.text(285, 0.2, '1', color=colors[1], fontsize=fs, ha='left', va='top')
axisD.text(285, 0.28, '2', color=colors[1], fontsize=fs, ha='left', va='top')
axisD.text(290, 0.39, '3', color=colors[1], fontsize=fs, ha='left', va='top')
axisD.text(315, 0.46, '4', color=colors[1], fontsize=fs, ha='left', va='top')
axisD.text(320, 0.6, '5', color=colors[1], fontsize=fs, ha='left', va='top')
aging_extras.draw_cell(axisD, 488, 521, 0.9, 16.7, y2xscale=1/333, toxic=True,
                                                           arrow=3.3, repair=2)
aging_extras.draw_const_env(axisD, 554, 0.8, 621, 1.0, y2xscale=1/333)
axisD.set_xlim([285, 630])
axisD.set_ylim([0.0, 1.0])
axisD.set_title('Age & Size Distributions')
left = 420
axisD.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
axisD.plot([left + 10], [asym_y], 'o', color=colorsD[0], markeredgecolor='none', markersize=ms)
axisD.text(left + 20, asym_y, 'Asymmetric, Optimal repair', va='center', ha='left', fontsize=fs)
axisD.plot([left + 10], [sym_y], 'o', color=colorsD[1], markeredgecolor='none', markersize=ms)
axisD.text(left + 20, sym_y, 'Symmetric, Optimal repair', va='center', ha='left', fontsize=fs)
setp( axisD.get_yticklabels(), visible=False)


axisCD.set_xlabel(r'Cellular total biomass ($P_{act} + P_{dam}$), (fg)', labelpad=22)

fig.process_subplots()
#axisB.tick_params(left='off')
#axisD.tick_params(left='off')
axisCD.tick_params(top="off", bottom="off", right="off", left="off")
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)
