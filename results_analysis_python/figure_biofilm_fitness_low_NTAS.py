#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results
import convert_fitness
import calc_roughness
import matplotlib.patches as mpatches

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'biofilms', 'preliminary', 'optimal_beta', 'Low')
output_path = os.path.join(input_path, 'low_NT_AS_fitness_strategy.pdf') #remember to change the end of the file name here!

fig = toolbox_plotting.BmcFigure(double_column=True) #, height='double')
# 5 cols instead - NR 006 007 008 DR
colors = ['red', '#009ACD', '#a203db', '#22c908', '#fbbf07']
# red, blue, purple, teal, yellow

ms = 1.25
fs = 6

#asymmetric, toxic
#test code with just these first, and then add the others in
axisASTP = fig.add_subplot('', 231) #plot population against time
axisASTM = fig.add_subplot('', 232) #plot biomass against time
axisASTG = fig.add_subplot('', 234) #plot growthrate against time
axisASTPG = fig.add_subplot('', 235) #plot growthrate/population against time
axislegend = fig.add_subplot('', 233, frameon = False) #legend
axisroughness = fig.add_subplot('', 236, frameon = False) #roughness

#pathsAST = ['no repair.xml', 'B=0.06.xml', 'B=0.07.xml', 'B=0.08.xml', 'DB.xml']

pathsAST = ['NTASANRL50_3day(20160222_1237)', 'NTASA006L50_3day(20160222_1203)', 'NTASA007L50_3day(20160222_1216)', 'NTASA008L50_3day(20160222_1216)', 'NTASADBL50_3day(20160222_1235)']
attributes = {'name':'population', 'name2':'mass', 'name3':'growthRate', 'header':'time,value,value2,value3'}
roughattributes = {'name':'Sigmaf', 'name2':'Sigma', 'name3':'Xf', 'name4':'Pf', 'header':'value,value2,value3,value4'}
allroughness = []
for i in range(5):
    path_input = os.path.join(input_path, pathsAST[i])
    convert_fitness.all_cells_growth_rate(path_input)
    results = os.path.join(input_path, pathsAST[i], 'fitness.xml')
    results_output = toolbox_results.ResultsOutput(path=results)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    pop = [float(r.vars['value']) for r in result_set.members]
    mass = [float(r.vars['value2']) for r in result_set.members]
    growth = [float(r.vars['value3']) for r in result_set.members]
    calc_roughness.calc_roughness(path_input)
    path_roughness = os.path.join(input_path, pathsAST[i], 'roughness.xml')
    roughness_output = toolbox_results.ResultsOutput(path=path_roughness)
    roughness_set = toolbox_results.ResultSet(roughness_output, roughattributes)
    roughness = [float(r.vars['value']) for r in roughness_set.members]
    allroughness.append(roughness[0])
    specgrowth = []
    for j in range(len(t)):
        specgrowth.append(pop[j]/mass[j])
    
    axisASTP.plot(t, pop, colors[i])
    axisASTM.plot(t, mass, colors[i])
    axisASTG.plot(t, growth, colors[i])
    axisASTPG.plot(t, specgrowth, colors[i])
    #setp( axisASTP.get_yticklabels(), visible=False)
    #setp( axisASTM.get_yticklabels(), visible=False)
    #setp( axisASTG.get_yticklabels(), visible=False)
    #setp( axisASTPG.get_yticklabels(), visible=False)
    axisASTP.set_xticks([0,24,48,72])
    axisASTM.set_xticks([0,24,48,72])
    axisASTG.set_xticks([0,24,48,72])
    axisASTPG.set_xticks([0,24,48,72])
    axisASTP.set_title('Population')
    axisASTM.set_title('Biomass')
    axisASTG.set_title('Growth Rate')
    axisASTPG.set_title('Specific Growth Rate')
no_repair = mpatches.Patch(color=colors[0], label='No repair')
b006 = mpatches.Patch(color=colors[1], label=r'$ \beta $= 0.06')
b007 = mpatches.Patch(color=colors[2], label=r'$ \beta $= 0.07')
b008 = mpatches.Patch(color=colors[3], label=r'$ \beta $= 0.08')
adaptive = mpatches.Patch(color=colors[4], label='Adaptive repair')
axislegend.legend(handles=[no_repair, b006, b007, b008, adaptive])
sigma = "$\sigma_f$ = "
no_repair_rough = mpatches.Patch(color=colors[0], label=(sigma +str(allroughness[0])))
b006_rough = mpatches.Patch(color=colors[1], label=(sigma +str(allroughness[1])))
b007_rough = mpatches.Patch(color=colors[2], label=(sigma +str(allroughness[2])))
b008_rough = mpatches.Patch(color=colors[3], label=(sigma +str(allroughness[3])))
adaptive_rough = mpatches.Patch(color=colors[4], label=(sigma +str(allroughness[4])))
axisroughness.legend(handles=[no_repair_rough, b006_rough, b007_rough, b008_rough, adaptive_rough])

#axisASTP.set_xticks([0, 24, 48, 72])
#axisAST.set_title('Asymmetric, Toxic')



fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.1, hspace=0.35)
fig.save(output_path)