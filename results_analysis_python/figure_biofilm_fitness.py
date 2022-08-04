#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import matplotlib.pyplot
import toolbox_basic
import toolbox_plotting
import toolbox_results
import convert_fitness
#import calc_roughness
import matplotlib.patches as mpatches

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'biofilms_used_for_optimal_comps', 'toxic', 'roughness_003556_symm')
output_path = os.path.join(input_path, 'T_M_S_fitness_fixed_beta_only.pdf') #remember to change the end of the file name here!

#_fixed_beta_only

fig = toolbox_plotting.BmcFigure(double_column=True) #, height='double')
# 5 cols instead - NR 006 007 008 DR
colors = ['red', '#009ACD', '#a203db', '#22c908', '#fbbf07']
# red, blue, purple, teal, yellow

ms = 1.25
fs = 6

#asymmetric, toxic
#test code with just these first, and then add the others in
axisASTP = fig.add_subplot('', 331) #plot population against time
axisASTM = fig.add_subplot('', 332) #plot biomass against time
axisASTG = fig.add_subplot('', 334) #plot growthrate against time
axisASTPG = fig.add_subplot('', 335) #plot growthrate/population against time
axistitle = fig.add_subplot('',342, frameon = False)
axislegend = fig.add_subplot('', 176, frameon = False) #legend

#pathsAST = ['B=0.06.xml', 'B=0.07.xml', 'B=0.08.xml', 'no repair.xml', 'DB.xml']

pathsAST = ['toxicitytrue_roughness0p003556_Symm0p0_betarepair0p06(20160303_1616)', 'toxicitytrue_roughness0p003556_Symm0p0_betarepair0p07(20160303_1705)', 'toxicitytrue_roughness0p003556_Symm0p0_betarepair0p08(20160303_1755)', 'toxicitytrue_roughness0p003556_Symm0p0_betarepair0p0(20160301_1021)', 'toxicitytrue_Symm0p0_roughness0p003556(20160301_0909)']
#pathsAST = ['toxicitytrue_Symm0p0_roughness0p003556(20160301_1726)']
attributes = {'name':'population', 'name2':'mass', 'name3':'growthRate', 'name4':'specificGrowth', 'header':'time,value,value2,value3,value4'}
for i in range(3):
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
    specgrowth = [float(r.vars['value4']) for r in result_set.members]
    for j in range(len(t)):
        mass[j] = mass[j]/1000000
        pop[j] = pop[j]/1000
        growth[j] = growth[j]/1000
    
    axisASTP.plot(t, pop, colors[i])
    axisASTP.set_ylabel('x1000', fontsize=6)
    axisASTM.plot(t, mass, colors[i])
    axisASTM.set_ylabel('ng', fontsize=6)
    axisASTG.plot(t, growth, colors[i])
    axisASTG.set_ylabel(r'pg $h^{-1}$', fontsize=6)
    axisASTPG.plot(t, specgrowth, colors[i])
    axisASTPG.set_ylabel(r'fg $h^{-1}$', fontsize = 6)
    setp(axisASTP.get_yticklabels(), fontsize=6)
    setp(axisASTM.get_yticklabels(), fontsize=6)
    setp(axisASTG.get_yticklabels(), fontsize=6)
    setp(axisASTPG.get_yticklabels(), fontsize=6)
    setp(axisASTP.get_xticklabels(), fontsize=6)
    setp(axisASTM.get_xticklabels(), fontsize=6)
    setp(axisASTG.get_xticklabels(), fontsize=6)
    setp(axisASTPG.get_xticklabels(), fontsize=6)
    axisASTP.set_xticks([0,24,48,72,96,120,144,168,192,216,240])
    axisASTP.set_xlabel('Time (h)', fontsize=7)
    axisASTM.set_xticks([0,24,48,72,96,120,144,168,192,216,240])
    axisASTM.set_xlabel('Time (h)', fontsize=7)
    axisASTG.set_xticks([0,24,48,72,96,120,144,168,192,216,240])
    axisASTG.set_xlabel('Time (h)', fontsize=7)
    axisASTPG.set_xticks([0,24,48,72,96,120,144,168,192,216,240])
    axisASTPG.set_xlabel('Time (h)', fontsize=7)
    axisASTP.set_title('Population', fontsize=7)
    axisASTM.set_title('Biomass', fontsize=7)
    axisASTG.set_title('Growth Rate', fontsize=7)
    axisASTPG.set_title('Specific Growth Rate', fontsize=7)
    #axistitle.set_title('Rough, Non-toxic, \n Symmetric', fontsize=8)
b006 = mpatches.Patch(color=colors[0], label=r'$ \beta $= 0.06')
b007 = mpatches.Patch(color=colors[1], label=r'$ \beta $= 0.07')
b008 = mpatches.Patch(color=colors[2], label=r'$ \beta $= 0.08')
no_repair = mpatches.Patch(color=colors[3], label='No repair')
adaptive = mpatches.Patch(color=colors[4], label='Adaptive repair')
axislegend.legend(handles=[b006, b007, b008, no_repair, adaptive], fontsize=7)

#axisASTP.set_xticks([0, 24, 48, 72])
#axisAST.set_title('Asymmetric, Toxic')


fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.4, hspace=0.4)
fig.save(output_path)