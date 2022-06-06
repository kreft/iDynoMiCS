#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results

#When using this plot, change the label_pos in toolbox_plotting to (0, 1.10)

#base_path = os.path.join('~', 'Dropbox', 'EclipseWorkspace', 'iDynoAge')
#base_path = os.path.join('C:\\', 'Users', 'RJW598', 'git', 'iDynoMiCS')
base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'comparison_clegg', 'Final_versions', 'Reproducing_figure_4')
output_path = os.path.join(input_path, 'figure_4_life_history_DEF.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True) #, height='double')
colors = ['red','#F781F3','blue','#81BEF7', '#fbbf07', 'green']

axisABC = fig.add_subplot('', 133, frameon=False)

ms = 1.25
fs = 6
axisA = fig.add_subplot('A', 332)
paths = ['T010AS12NR100_agent_State_last.xml', 'T010S12NR100_agent_State_last.xml']
         #, 'sn10_const_toxic_agent_state_last.xml']
colorsA = [colors[0], colors[2]]
for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    agent_output = toolbox_results.AgentOutput(path=results_path)
    species = toolbox_results.SpeciesOutput(agent_output)
    aging_extras.scatter_population(axisA, species, 'totalBiomass', 'age',
                                                color=colorsA[i], markersize=ms)
width, height = 5, 0.12
left, right, bottom = 285, 627, -0.01
'''
axisA.plot([left + width, left, left, left + width],
                                          [bottom]*2 + [bottom+height]*2, 'k-')
axisA.plot([right - width, right, right, right - width],
                                          [bottom]*2 + [bottom+height]*2, 'k-')
'''
axisA.text(290, 0.10, 'Generation 0',
                            color=colors[0], fontsize=fs, ha='left', va='top')
axisA.text(290, 0.27, '1', color=colors[0], fontsize=fs, ha='left', va='top')
axisA.text(290, 0.45, '2', color=colors[0], fontsize=fs, ha='left', va='top')
axisA.text(290, 0.76, '3', color=colors[0], fontsize=fs, ha='left', va='top')
axisA.text(290, 0.97, '4', color=colors[0], fontsize=fs, ha='left', va='top')
aging_extras.draw_cell(axisA, 488, 521, 0.9, 16.7, y2xscale=1/333, toxic=True,
                                                           arrow=3.3, repair=0)
aging_extras.draw_const_env(axisA, 554, 0.8, 621, 1.0, y2xscale=1/333)
axisA.set_xlim([285, 630])
axisA.set_ylim([0.0, 1.0])
axisA.set_xticks([300, 400, 500, 600])
#axisA.set_title('Age & Size Distributions')
axisA.set_ylabel(r'Cellular age $P_{dam}/(P_{act}+P_{dam})$', fontsize=8)
# Legend
left, right = 451, 630
top, bottom = 0.78, 0.62
asym_y = 0.735
sym_y = 0.665
axisA.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
axisA.plot([left + 10], [asym_y], 'o', color=colorsA[0], markeredgecolor='none', markersize=3)
axisA.text(left + 20, asym_y, 'Asymmetric, No repair', va='center', ha='left', fontsize=fs)
axisA.plot([left + 10], [sym_y], 'o', color=colorsA[1], markeredgecolor='none', markersize=3)
axisA.text(left + 20, sym_y, 'Symmetric, No repair', va='center', ha='left', fontsize=fs)


axisB = fig.add_subplot('B', 335)
paths = ['T010AS12007100_agent_State_last.xml', 'T010S12007100_agent_State_last.xml']
         #, 'so10_const_toxic_agent_state_last.xml']
colorsB = [colors[1], colors[3]]
for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    agent_output = toolbox_results.AgentOutput(path=results_path)
    species = toolbox_results.SpeciesOutput(agent_output)
    aging_extras.scatter_population(axisB, species, 'totalBiomass', 'age', color=colorsB[i], markersize=ms)
axisB.text(285, 0.05, '0', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(285, 0.2, '1', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(285, 0.33, '2', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(290, 0.42, '3', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(305, 0.49, '4', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(320, 0.6, '5', color=colors[1], fontsize=fs, ha='left', va='top')
axisB.text(340, 0.7, '6', color=colors[1], fontsize=fs, ha='left', va='top')
aging_extras.draw_cell(axisB, 488, 521, 0.9, 16.7, y2xscale=1/333, toxic=True,
                                                           arrow=3.3, repair=2)
aging_extras.draw_const_env(axisB, 554, 0.8, 621, 1.0, y2xscale=1/333)
axisB.set_xlim([285, 630])
axisB.set_ylim([0.0, 1.0])
axisB.set_xticks([300, 400, 500, 600])
axisB.set_title('Age & Size Distributions', fontsize=9)
axisB.set_xlabel('Cellular total biomass (fg)', fontsize=8)
left = 420
right = 630
axisB.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
axisB.plot([left + 10], [asym_y], 'o', color=colorsB[0], markeredgecolor='none', markersize=3)
axisB.text(left + 20, asym_y, 'Asymmetric, Optimal repair', va='center', ha='left', fontsize=fs)
axisB.plot([left + 10], [sym_y], 'o', color=colorsB[1], markeredgecolor='none', markersize=3)
axisB.text(left + 20, sym_y, 'Symmetric, Optimal repair', va='center', ha='left', fontsize=fs)
setp( axisB.get_xticklabels(), visible=False)
axisB.set_ylabel(r'Cellular age $P_{dam}/(P_{act}+P_{dam})$', fontsize=8)



axisC = fig.add_subplot('C', 338)
paths = ['T010AS12OB100_agent_State_last.xml', 'T010S12OB100_agent_State_last.xml']
         #, '.xml']
colorsC = [colors[4], colors[5]]
for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    agent_output = toolbox_results.AgentOutput(path=results_path)
    species = toolbox_results.SpeciesOutput(agent_output)
    aging_extras.scatter_population(axisC, species, 'totalBiomass', 'age', color=colorsC[i], markersize=ms)
axisC.text(285, 0.05, '0', color=colors[4], fontsize=fs, ha='left', va='top')
axisC.text(280, 0.2, '1', color=colors[4], fontsize=fs, ha='left', va='top')
axisC.text(285, 0.35, '2', color=colors[4], fontsize=fs, ha='left', va='top')
axisC.text(300, 0.50, '3+', color=colors[4], fontsize=fs, ha='left', va='top')
aging_extras.draw_cell(axisC, 488, 521, 0.9, 16.7, y2xscale=1/333, toxic=True,
                                                           arrow=3.3, repair=2)
aging_extras.draw_const_env(axisC, 554, 0.8, 621, 1.0, y2xscale=1/333)
axisC.set_xlim([285, 630])
axisC.set_ylim([0.0, 1.0])
axisC.set_xticks([300, 400, 500, 600])
#axisC.set_title('Age & Size Distributions')
left = 420
right = 630
axisC.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
axisC.plot([left + 10], [asym_y], 'o', color=colorsC[0], markeredgecolor='none', markersize=3)
axisC.text(left + 20, asym_y, 'Asymmetric, Adaptive repair', va='center', ha='left', fontsize=fs)
axisC.plot([left + 10], [sym_y], 'o', color=colorsC[1], markeredgecolor='none', markersize=3)
axisC.text(left + 20, sym_y, 'Symmetric, Adaptive repair', va='center', ha='left', fontsize=fs)
setp( axisC.get_xticklabels(), visible=False)
axisC.set_ylabel(r'Cellular age $P_{dam}/(P_{act}+P_{dam})$', fontsize=8)

#move this to the bottom if you actually want there to be a sub plot D!
fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.0)
fig.save(output_path)