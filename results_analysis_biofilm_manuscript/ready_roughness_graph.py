#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import toolbox_basic
import toolbox_results
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

inp_path = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/roughness/roughness_figure/'

fig = plt.figure(figsize=(8.27, 4))
colors = ['r','b','g']

axisABC = plt.subplot(111)
high = ['NANRH152_roughness_1.xml', 'NANRH152_roughness_2.xml', 'NANRH152_roughness_3.xml', 'NANRH154_roughness_4.xml']
medium = ['NANRM154_roughness_1.xml', 'NANRM152_roughness_2.xml', 'NANRM152_roughness_3.xml', 'NANRM152_roughness_4.xml']
low = ['NANRL152_roughness_1.xml', 'NANRL152_roughness_2.xml', 'NANRL152_roughness_3.xml', 'NANRL154_roughness_4.xml']
attributes = {'name':'Sigmaf', 'name2':'Sigma', 'name3':'Xf', 'name4':'Pf', 'name5':'height', 'header':'value,value2,value3,value4,value5'}

alpha = 0.2
for a in range(3):
    if a == 0: paths = high
    if a == 1: paths = medium
    if a == 2: paths = low
    maxi = 0
    sortbyheight = numpy.zeros((65, 200), dtype=float)
    countrow = [0]*65
    mini = 100
    for b in range(4): 
        sim = os.path.join(inp_path, paths[b])
        results_output = toolbox_results.ResultsOutput(path=sim)
        result_set = toolbox_results.ResultSet(results_output, attributes)
        if b == 0:
            rough0 = [float(r.vars['value']) for r in result_set.members]
            high0 = [float(r.vars['value5']) for r in result_set.members]
            rough, high = rough0, high0
        if b == 1:
            rough1 = [float(r.vars['value']) for r in result_set.members]
            high1 = [float(r.vars['value5']) for r in result_set.members]
            rough, high = rough1, high1
        if b == 2:
            rough2 = [float(r.vars['value']) for r in result_set.members]
            high2 = [float(r.vars['value5']) for r in result_set.members]
            rough, high = rough2, high2
        if b == 3:
            rough3 = [float(r.vars['value']) for r in result_set.members]
            high3 = [float(r.vars['value5']) for r in result_set.members]
            rough, high = rough3, high3
        for c in range(len(high)):
            high[c] = int(high[c])
            row = high[c]
            column = countrow[high[c]]
            countrow[high[c]]+=1
            if not rough[c] == 0: 
                sortbyheight[row, column] = rough[c]
    roughness = []
    standev = []
    height = []
    for d in range(65):
        numberstouse = []
        for e in range(20):
            if d > 20 and sortbyheight[d,e]==0:
                continue
            sortbyheight[d,e]
            numberstouse.append(sortbyheight[d,e])
        if len(numberstouse) == 1:
            roughness.append(numberstouse[0])
            standev.append(0)
            height.append(d*4)
        if len(numberstouse)>0 and not len(numberstouse) == 1:
            roughness.append(numpy.mean(numberstouse))
            standev.append(numpy.std(numberstouse))
            height.append(d*4)
    stdev1 = []
    stdev2 = []
    heightlow = []
    heightup = []
    for f in range(len(standev)):
        stdev1.append(roughness[f]-standev[f])
        stdev2.append(roughness[f]+standev[f])
        heightlow.append(height[f]-2)
        heightup.append(height[f]+2)
    axisABC.fill_betweenx(height, stdev1, stdev2, facecolor = colors[a], alpha=alpha)
    labelhigh = (r'$\delta^2$ = 0.014222')
    labelmedium = (r'$\delta^2$ = 0.003556')
    labellow = (r'$\delta^2$ = 0.000889')
    line = axisABC.plot(roughness, height, color=colors[a])
patch1 = mpatches.Patch(color=colors[0], label=labelhigh)
patch2 = mpatches.Patch(color=colors[1], label=labelmedium)
patch3 = mpatches.Patch(color=colors[2], label=labellow)
plt.legend(handles=[patch1, patch2, patch3], fontsize=9, loc='upper right')
 
os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/roughness/')      
axisABC.set_xticks([0, 1, 2, 3, 4, 5, 5, 6])
axisABC.set_ylim([0, 180])
axisABC.set_xlim([0, 6])
axisABC.set_ylabel(r'Biofilm height ($\mu$m)')
axisABC.set_xlabel(r'Roughness ($\sigma_{f}$)')
fig.savefig('Roughness_graph.png', bbox_inches='tight', dpi=600)