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
import convert_death
import convert_fitness
#import calc_roughness
import matplotlib.patches as mpatches

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'feb_2017')
output_path = os.path.join(input_path, 'death_test_biom.pdf') #remember to change the end of the file name here!

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

pathsAST = ['bio_comp_16_cells_1_biom(20170205_1725)']
#pathsAST = ['toxicitytrue_Symm0p0_roughness0p003556(20160301_1726)']
attributes = {'name':'population', 'name2':'mass', 'name3':'growthRate', 'name4':'specificGrowth', 'header':'time,value,value2,value3,value4'}
for i in range(1):
    OldieAdetachmenttime, OldieAdetachmentpop, OldieAdetachmentmass, OldieAdetachmentgrowth, OldieAdetachmentspecgrowth = [], [], [], [], []
    OldieBdetachmenttime, OldieBdetachmentpop, OldieBdetachmentmass, OldieBdetachmentgrowth, OldieBdetachmentspecgrowth = [], [], [], [], []
    OldieAtime, OldieApop, OldieAmass, OldieAgrowth, OldieAspecgrowth = [], [], [], [], []
    OldieBtime, OldieBpop, OldieBmass, OldieBgrowth, OldieBspecgrowth = [], [], [], [], []
    times = [OldieAdetachmenttime, OldieBdetachmenttime, OldieAtime, OldieBtime]
    pops = [OldieAdetachmentpop, OldieBdetachmentpop, OldieApop, OldieBpop]
    masses = [OldieAdetachmentmass, OldieBdetachmentmass, OldieAmass, OldieBmass]
    growths = [OldieAdetachmentgrowth, OldieBdetachmentgrowth, OldieAgrowth, OldieBgrowth]
    specgrowths = [OldieAdetachmentspecgrowth, OldieBdetachmentspecgrowth, OldieAspecgrowth, OldieBspecgrowth]
    path_input = os.path.join(input_path, pathsAST[i])
    convert_death.all_cells_death_rate(path_input)
    convert_fitness.all_cells_growth_rate(path_input)
    for j in range(4):
        if j == 0:
            resultsusing = 'OldieAdetachment.xml'
        elif j == 1:
            resultsusing = 'OldieBdetachment.xml'
        elif j == 2:
            resultsusing = 'OldieAfitness.xml'
        elif j == 3:
            resultsusing = 'OldieBfitness.xml'
        results = os.path.join(input_path, pathsAST[i], resultsusing)
        results_output = toolbox_results.ResultsOutput(path=results)
        result_set = toolbox_results.ResultSet(results_output, attributes)
        for result in result_set.members:
            result.vars['time'] = float(result.vars['time'])
        result_set.members.sort(key=lambda r: r.vars['time'])
        for r in result_set.members:
            t = float(r.vars['time'])
            times[j].append(t)
            pop = float(r.vars['value'])
            pops[j].append(pop)
            mass = float(r.vars['value2'])
            masses[j].append(mass)
            growth = float(r.vars['value3'])
            growths[j].append(growth)
            specgrowth = float(r.vars['value4'])
            specgrowths[j].append(specgrowth)
    for a in range(len(pops[0])):
        if a > 0:
            pops[0][a] += pops[0][a-1]
            pops[1][a] += pops[1][a-1]
            masses[0][a] += masses[0][a-1]
            masses[1][a] += masses[1][a-1]
        
    
    for k in range(2): 
        if k == 0:
            t, pop, mass, t1, detachpop, detachmass = OldieAtime, OldieApop, OldieAmass, OldieAdetachmenttime, OldieAdetachmentpop, OldieAdetachmentmass
        elif k == 1:
            t, pop, mass, t1, detachpop, detachmass = OldieBtime, OldieBpop, OldieBmass, OldieBdetachmenttime, OldieBdetachmentpop, OldieBdetachmentmass
        print len(t)
        print len(pop)
        axisASTP.plot(t, pop, colors[k]) #pop vs time
        axisASTP.set_ylabel('', fontsize=6)
        axisASTM.plot(t, mass, colors[k]) #biomass vs time
        axisASTM.set_ylabel('', fontsize=6)
        axisASTG.plot(t1, detachpop, colors[k]) #pop detachment vs time
        axisASTG.set_ylabel(r'', fontsize=6)
        axisASTPG.plot(t1, detachmass, colors[k]) #detachment biomass vs time
        axisASTPG.set_ylabel(r'', fontsize = 6)
        setp(axisASTP.get_yticklabels(), fontsize=6)
        setp(axisASTM.get_yticklabels(), fontsize=6)
        setp(axisASTG.get_yticklabels(), fontsize=6)
        setp(axisASTPG.get_yticklabels(), fontsize=6)
        setp(axisASTP.get_xticklabels(), fontsize=6)
        setp(axisASTM.get_xticklabels(), fontsize=6)
        setp(axisASTG.get_xticklabels(), fontsize=6)
        setp(axisASTPG.get_xticklabels(), fontsize=6)
        axisASTP.set_xlabel('Time (h)', fontsize=7)
        axisASTM.set_xlabel('Time (h)', fontsize=7)
        axisASTG.set_xlabel('Time (h)', fontsize=7)
        axisASTPG.set_xlabel('Time (h)', fontsize=7)
        axisASTP.set_title('Population', fontsize=7)
        axisASTM.set_title('Biomass', fontsize=7)
        axisASTG.set_title('Detachment', fontsize=7)
        axisASTPG.set_title('Detached biomass', fontsize=7)
    #axistitle.set_title('Rough, Non-toxic, \n Symmetric', fontsize=8)
b006 = mpatches.Patch(color=colors[0], label='Adaptive repair')
b007 = mpatches.Patch(color=colors[1], label='Damage segregation')
axislegend.legend(handles=[b006, b007], fontsize=7)

#axisASTP.set_xticks([0, 24, 48, 72])
#axisAST.set_title('Asymmetric, Toxic')


fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.4, hspace=0.4)
fig.save(output_path)