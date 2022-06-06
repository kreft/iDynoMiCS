#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'Figure_4')
output_path = os.path.join(input_path, 'figure_4_life_history_Prep_Ptot_sameaxis.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')
ms = 2

axisA = fig.add_subplot('B', 221)
paths = ['T010AS12OBConstEnv1_0D_Ind_results_for_Prep_Ptot.xml']
colors = ['#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#06e4ea', '#0066FF', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea']
attributes = {'name':'age', 'name2':'activeBiomassGrowth', 'name3':'activeBiomassRepair', 'name4':'inactiveBiomassGrowth',
'name5':'inactiveBiomassRepair', 'name6':'generation', 'header':'time,value,value2,value3,value4,value5,value6'}

muS = 0.6
Yr = 0.8

for i in range(1):    
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    age = [float(r.vars['value']) for r in result_set.members]
    actGro = [float(r.vars['value2']) for r in result_set.members]
    actRep = [float(r.vars['value3']) for r in result_set.members]
    inactGro = [float(r.vars['value4']) for r in result_set.members]
    inactRep = [float(r.vars['value5']) for r in result_set.members]
    gen = [int(r.vars['value6']) for r in result_set.members]
    linex = ([0,0.5])
    liney = ([0,0.5])
    x = []
    y = []
    oldx = []
    oldy = []
    color = []
    axisA.plot(linex, liney, 'k')
    current_gen = 0
    count = 0
    for j in range(0, len(age)):
            pActGro = float(actGro[j])
            pActRep = float(actRep[j])
            pInactGro = float(inactGro[j])
            pInactRep = float(inactRep[j])
            pRep = (pActRep + pInactRep)
            pTot = (pActGro + pActRep + pInactGro + pInactRep)
            if gen[j] > current_gen:
                axisA.plot(x, y, colors[current_gen])
                oldx.append(x[-1])
                oldy.append(y[-1])
                x = []
                y = []
            if current_gen < 7:
                x.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS) * (1 / (1 - age[j]))) - 1))])
                y.append([pRep/pTot])
                if gen[j] > current_gen:
                    w = x[0]
                    wx = ([w, oldx[0]])
                    z = y[0]
                    zy = ([z, oldy[0]])
                    axisA.plot(wx, zy, linestyle='--', color='#565353' )
                    axisA.plot(oldx[0], oldy[0], marker='o', color=colors[current_gen])
                    current_gen += 1
                    oldx = []
                    oldy = []
            else:
                break
fs = 8 
axisA.set_xlabel(r'$ \^\beta $', style='italic')
axisA.set_ylabel(r'$P_{rep}/P_{tot}$')
#axisA.set_title('Following old pole cell over time (Asymmetric division)')
axisA.set_xlim([0, 0.45])
axisA.set_ylim([0, 0.45])
#axisA.text(0.175,0.165,'x=y', color='k', fontsize=10, style='italic')


axisB = fig.add_subplot('D', 223)
paths = ['T010S12OBConstEnvCleggComp_results_for_Prep_Ptot.xml']
attributes = {'name':'age', 'name2':'activeBiomassGrowth', 'name3':'activeBiomassRepair', 'name4':'inactiveBiomassGrowth',
'name5':'inactiveBiomassRepair', 'name6':'generation', 'header':'time,value,value2,value3,value4,value5,value6'}

muS = 0.6
Yr = 0.8

for i in range(1):    
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    age = [float(r.vars['value']) for r in result_set.members]
    actGro = [float(r.vars['value2']) for r in result_set.members]
    actRep = [float(r.vars['value3']) for r in result_set.members]
    inactGro = [float(r.vars['value4']) for r in result_set.members]
    inactRep = [float(r.vars['value5']) for r in result_set.members]
    gen = [int(r.vars['value6']) for r in result_set.members]
    linex = ([0,0.45])
    liney = ([0,0.45])
    x = []
    y = []
    oldx = []
    oldy = []
    color = []
    axisB.plot(linex, liney, 'k')
    current_gen = 0
    count = 0
    for j in range(0, len(age)):
            pActGro = float(actGro[j])
            pActRep = float(actRep[j])
            pInactGro = float(inactGro[j])
            pInactRep = float(inactRep[j])
            pRep = (pActRep + pInactRep)
            pTot = (pActGro + pActRep + pInactGro + pInactRep)
            if gen[j] > current_gen:
                axisB.plot(x, y, colors[current_gen])
                oldx.append(x[-1])
                oldy.append(y[-1])
                x = []
                y = []
            x.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS) * (1 / (1 - age[j]))) - 1))])
            y.append([pRep/pTot])
            if gen[j] > current_gen:
                w = x[0]
                wx = ([w, oldx[0]])
                z = y[0]
                zy = ([z, oldy[0]])
                axisB.plot(w, z, marker='o', color=colors[current_gen])
                current_gen += 1
                oldx = []
                oldy = []
fs = 8 
axisB.set_xlabel(r'$ \^\beta $', style='italic')
axisB.set_ylabel(r'$P_{rep}/P_{tot}$')
#axisB.set_title('Following old pole cell over time (Symmetric division)')
axisB.set_xlim([0, 0.45])
axisB.set_ylim([0, 0.45])
#axisB.text(0.047,0.045,'x=y', color='k', fontsize=10, style='italic')


fig.process_subplots()
fig.subplots_adjust(left=0.11,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)