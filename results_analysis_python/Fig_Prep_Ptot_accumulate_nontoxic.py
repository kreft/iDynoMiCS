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
input_path = os.path.join(base_path, 'results', 'comparison_clegg', 'Final_versions', 'Reproducing_figure_4', 'nontoxic')
output_path = os.path.join(input_path, 'figure_4_life_history_Prep_Ptot_accumulate.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')
ms = 2

axisA = fig.add_subplot('A', 221)
paths = ['NT010AS12OB3day_results_for_Prep_Ptot.xml']
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
        if t[j] < 25:
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
            if current_gen <= 7:
                x.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS)) - 1))])
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
axisA.set_title('Asymmetric')


axisB = fig.add_subplot('B', 222)
paths = ['NT010S12OB3day_results_for_Prep_Ptot.xml']
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
    linew = ([0,0.5])
    linez = ([0,0.5])
    w = []
    z = []
    oldw = []
    oldz = []
    color = []
    axisB.plot(linew, linez, 'k')
    current_gen = 0
    count = 0
    for j in range(0, len(age)):
            if t[j] < 25:
                pActGro = float(actGro[j])
                pActRep = float(actRep[j])
                pInactGro = float(inactGro[j])
                pInactRep = float(inactRep[j])
                pRep = (pActRep + pInactRep)
                pTot = (pActGro + pActRep + pInactGro + pInactRep)
                if gen[j] > current_gen:
                    axisB.plot(w, z, colors[current_gen])
                    oldw.append(w[-1])
                    oldz.append(z[-1])
                    w = []
                    z = []
                if current_gen <= 7:
                    w.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS)) - 1))])
                    z.append([pRep/pTot])
                    if gen[j] > current_gen:
                        m = w[0]
                        mw = ([m, oldw[0]])
                        n = z[0]
                        nz = ([n, oldz[0]])
                        axisB.plot(mw, nz, linestyle='--', color='#565353')
                        axisB.plot(oldw[0], oldz[0], marker='o', color=colors[current_gen])
                        current_gen += 1
                        oldw = []
                        oldz = []
                else:
                    break
fs = 8 
axisB.set_xlabel(r'$ \^\beta $', style='italic')
axisB.set_xlim([0, 0.055])
axisB.set_ylim([0, 0.055])
axisB.set_title('Symmetric')

axisC = fig.add_subplot('C', 223)
paths = ['NT010AS12OB3day_results_for_Prep_Ptot.xml']
colors = ['#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#06e4ea', '#0066FF', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea',
          '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#06e4ea', '#0066FF', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea',
          '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#06e4ea', '#0066FF', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea']
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
    x = []
    y = []
    xvec = []
    yvec = []
    pReppTot = []
    beta = []
    for j in range(0, len(age)):
        if t[j] < 25:
            pActGro = float(actGro[j])
            pActRep = float(actRep[j])
            pInactGro = float(inactGro[j])
            pInactRep = float(inactRep[j])
            pRep = (pActRep + pInactRep)
            pTot = (pActGro + pActRep + pInactGro + pInactRep)
            pReppTot.append([pRep/pTot])
            beta.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS)) - 1))])
    current_gen = 0
    for k in range(0, len(gen)):
        if t[k] < 25:
            if gen[k] == current_gen:
                xvec.append(beta[k])
                yvec.append(pReppTot[k])
            elif gen[k] > current_gen:
                xi = sum(xvec)/len(xvec)
                yi = sum(yvec)/len(yvec)
                x.append(sum(xvec)/len(xvec))
                y.append(sum(yvec)/len(yvec))
                axisC.plot(xi, yi, 'o', color=colors[current_gen])
                current_gen += 1
        
    axisC.plot(x, y, '--', color='#565353')
    
fs = 8 
axisC.set_xlabel(r'$ \^\beta $', style='italic')
axisC.set_ylabel(r'$P_{rep}/P_{tot}$')
axisC.set_title('Asymmetric')
axisC.set_ylim([0, 0.15])

axisD = fig.add_subplot('D', 224)
paths = ['NT010S12OB3day_results_for_Prep_Ptot.xml']
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
    x = []
    y = []
    xvec = []
    yvec = []
    pReppTot = []
    beta = []
    for j in range(0, len(age)):
        if t[j] < 25:
            pActGro = float(actGro[j])
            pActRep = float(actRep[j])
            pInactGro = float(inactGro[j])
            pInactRep = float(inactRep[j])
            pRep = (pActRep + pInactRep)
            pTot = (pActGro + pActRep + pInactGro + pInactRep)
            pReppTot.append([pRep/pTot])
            beta.append([float(( age[j] / (1 - age[j]))  * ( math.sqrt((Yr / muS)) - 1))])
    current_gen = 0
    for k in range(0, len(gen)):
        if t[k] < 20:
            if gen[k] == current_gen:
                xvec.append(beta[k])
                yvec.append(pReppTot[k])
            elif gen[k] > current_gen:
                xi = sum(xvec)/len(xvec)
                yi = sum(yvec)/len(yvec)
                x.append(sum(xvec)/len(xvec))
                y.append(sum(yvec)/len(yvec))
                axisD.plot(xi, yi, 'o', color=colors[current_gen])
                current_gen += 1
        
    axisD.plot(x, y, '--', color='#565353')
    
fs = 8 
axisD.set_xlabel(r'$ \^\beta $', style='italic')
axisB.set_xlim([0, 0.055])
axisB.set_ylim([0, 0.055])
axisD.set_xticks([0, 0.005, 0.010, 0.015, 0.020, 0.025])
axisD.set_title('Symmetric')
#axisC.set_ylabel(r'$P_{rep}/P_{tot}$')





fig.process_subplots()
fig.subplots_adjust(left=0.11,right=0.99,top=0.96,bottom=0.08,wspace=0.25, hspace=0.35)
fig.save(output_path)