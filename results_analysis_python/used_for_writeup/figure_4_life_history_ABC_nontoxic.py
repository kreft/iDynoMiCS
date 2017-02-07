#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results
from mpl_toolkits.axes_grid1 import make_axes_locatable

#base_path = os.path.join('~', 'Dropbox', 'EclipseWorkspace', 'iDynoAge')
#base_path = os.path.join('C:\\', 'Users', 'RJW598', 'git', 'iDynoMiCS')
base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'comparison_clegg', 'Final_versions', 'Reproducing_figure_4', 'nontoxic')
output_path = os.path.join(input_path, 'figure_4_life_history_ABC_3D_adapt.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')


axisA = fig.add_subplot('A', 221)
#SOB, SNR, S007, ASNR, AS007, ASOB
paths = ['NT010S12OB_life_history_results.xml', 'NT010S12NR_life_history_results.xml', 'NT010S120023_life_history_results.xml','NT010AS12NR_life_history_results.xml','NT010AS120023_life_history_results.xml', 'NT010AS12OB_life_history_results.xml']
colorsA = ['green', 'blue', '#81BEF7', 'red','#F781F3', '#fbbf07']
#attributes1 = {'name':'specific growth rate', 'header':'time,value'}
attributes = {'name':'specific growth rate', 'name2':'age', 'name3':'generation', 'header':'time,value,value2,value3'}

xlim, ylim = 26, 0.6

for i in range(6):    
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    val = [float(r.vars['value']) for r in result_set.members]
    gen = [float(r.vars['value3']) for r in result_set.members]
    j = 1
    while j < len(t)-1:
        x=[]
        y=[]
        while j < len(t)-1 and [gen[j]] == [gen[j+1]]:
            x.append([t[j+1]])
            y.append([val[j+1]])
            j += 1
        axisA.plot(x, y, color=colorsA[i])
        #axisA.plot(divx, divy, 'k*')
        j += 1
fs = 8 
axisA.text(1, 0.57, '0', color='#8B4513', fontsize=fs)
axisA.text(0.5, 0.47, '1', color='#8B4513', fontsize=fs)
axisA.text(4.25, 0.40,'2', color='#8B4513', fontsize=fs)
axisA.text(6, 0.33, '3', color='#8B4513', fontsize=fs)
axisA.text(7.5, 0.19, '4', color=colorsA[3], fontsize=fs)
axisA.text(11.5, 0.19, '5', color=colorsA[3], fontsize=fs)
axisA.text(17, 0.09, '6', color=colorsA[3], fontsize=fs)
axisA.text(21, 0.01, '7', color=colorsA[3], fontsize=fs)
axisA.text(7.5, 0.30, '4', color='#8B4513', fontsize=fs)
axisA.text(10.5, 0.29, '5', color='#8B4513', fontsize=fs)
axisA.text(13, 0.29, '6', color='#8B4513', fontsize=fs)
axisA.text(15.5, 0.29, '7', color='#8B4513', fontsize=fs)
axisA.text(17.5, 0.20, '8', color=colorsA[4], fontsize=fs)
axisA.text(21, 0.18, '9', color=colorsA[4], fontsize=fs)
axisA.text(23, 0.04, '10', color=colorsA[4], fontsize=fs)
axisA.text(17.25, 0.27, '8', color=colorsA[5], fontsize=fs)
axisA.text(20, 0.27, '9', color=colorsA[5], fontsize=fs)
axisA.text(22, 0.27, '10', color=colorsA[5], fontsize=fs)
axisA.text(23.5, 0.20, '11', color=colorsA[5], fontsize=fs)

aging_extras.draw_cell(axisA, 0.59*xlim, 0.69*xlim, 0.9*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=False, arrow=0.01*xlim)
aging_extras.draw_const_env(axisA, 0.79*xlim, 0.8*ylim, 0.99*xlim, ylim, y2xscale=ylim/xlim)
axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 4, 8, 12, 16, 20, 24])
#axisA.set_ylim([0, ylim])
axisA.set_xlabel(r'Time (h)')
axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
axisA.set_title('Following Old Pole Cell Over Time')


axisB = fig.add_subplot('B', 222)
paths = ['NT010AS12NR100_results.xml', 'NT010AS12007100_results.xml', 'NT010AS12OB100_results.xml']
attributes = {'name':'specific growth rate population structure',
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
    axisB.plot(x, y, color=colorsB[i])
    
axisB.text(0.0775, 0.57, '0', color=colorsB[0], fontsize=fs)
axisB.text(0.12, 0.51, '1', color=colorsB[0], fontsize=fs)
axisB.text(0.025, 0.40, '2-3', color=colorsB[0], fontsize=fs)
axisB.text(0.095, 0.52, '0', color=colorsB[1], fontsize=fs)
axisB.text(0.13, 0.47, '1', color=colorsB[1], fontsize=fs)
axisB.text(0.04, 0.425, '2-3', color=colorsB[1], fontsize=fs)
axisB.text(0.019, 0.375, '4-7', color=colorsB[1], fontsize=fs)
axisB.text(0.0815, 0.54, '0', color=colorsB[2], fontsize=fs)
axisB.text(0.159, 0.485, '1', color=colorsB[2], fontsize=fs)
axisB.text(0.03, 0.445, '2-3', color=colorsB[2], fontsize=fs)
axisB.text(0.015, 0.32, '4-', color=colorsB[2], fontsize=fs)

axisB.set_xlim([0, 0.17])
axisB.set_xticks([0, 0.04, 0.08, 0.12, 0.16])
axisB.set_ylim([0, 0.6])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate: Distribution')
setp( axisB.get_yticklabels(), visible=False)

pad = 0.00
side = "right"
divider = make_axes_locatable(axisB)
caxB = divider.append_axes(side, size="25%", pad=pad, axisbg='none', frameon=True)
caxB.set_xticklabels(['']*10)
caxB.set_yticklabels(['']*10)
caxB.spines['left'].set_color('none')
caxB.tick_params(top="off", right="off", bottom="off", left="off")
caxB.scatter(0.5, 0.495773, color=colorsB[0])
plt = ([0.5,0.5])
NR = ([0.495773+0.093651, 0.495773-0.093651])
caxB.plot(plt, NR, '--', color=colorsB[0])
caxB.plot([0.45, 0.55], [NR[0], NR[0]], color=colorsB[0])
caxB.plot([0.45, 0.55], [NR[1], NR[1]], color=colorsB[0])
caxB.scatter(1, 0.484955, color=colorsB[1])
plt = ([1,1])
opt = ([0.484955+0.055283, 0.484955-0.055283])
caxB.plot(plt, opt, '--', color=colorsB[1])
caxB.plot([0.95, 1.05], [opt[0], opt[0]], color=colorsB[1])
caxB.plot([0.95, 1.05], [opt[1], opt[1]], color=colorsB[1])
caxB.scatter(1.5, 0.497153, color=colorsB[2])
plt = ([1.5, 1.5])
AR = ([0.497153+0.076518, 0.497153-0.076518])
caxB.plot(plt, AR, '--', color=colorsB[2])
caxB.plot([1.45, 1.55], [AR[0], AR[0]], color=colorsB[2])
caxB.plot([1.45, 1.55], [AR[1], AR[1]], color=colorsB[2])
caxB.set_xlim([0, 2])
caxB.set_ylim([0, 0.6])
caxB.set_title(r'Means $\pm$ SD')


axisC = fig.add_subplot('C', 212)
paths = ['NT010AS120023_life_history_results.xml', 'NT010S120023_life_history_results.xml']
pathsprot = ['NT010AS12OB3day_results_for_Prep_Ptot.xml', 'NT010S12OB3day_results_for_Prep_Ptot.xml']
colorsC = ['#fbbf07', 'green', '#81BEF7']
attributes = {'name':'specific growth rate', 'name2':'age', 'name3':'generation', 'header':'time,value,value2,value3'}
attributesprot = {'name':'age', 'name2':'activeBiomassGrowth', 'name3':'activeBiomassRepair', 'name4':'inactiveBiomassGrowth',
'name5':'inactiveBiomassRepair', 'name6':'generation', 'header':'time,value,value2,value3,value4,value5,value6'}

muS = 0.6
Yr = 0.8

for i in range(2):
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    #updated so as this takes value 1 - for my results this is specific growth
    #rate, which value 2 is age
    val = [float(r.vars['value2']) for r in result_set.members]
    gen = [float(r.vars['value3']) for r in result_set.members]
    growth = [float(r.vars['value']) for r in result_set.members]
    j = 0
    while j < len(t)-1:
        x, y = [], []
        while (j < len(t)-1) and [gen[j]] == [gen[j+1]]:    
            age = val[j+1]
            x.append(t[j+1])
            y.append([float(( age / (1 - age))  * ( math.sqrt((Yr / muS)) - 1))])
            j += 1
        axisC.plot(x, y, color=colorsC[i])
        j += 1       
    fs = 8 
xall = [0, 72]
yall = [0.02, 0.02]
axisC.plot(xall, yall, color=colorsC[2])
for i in range(2):
    results_path = os.path.join(input_path, pathsprot[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributesprot)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    Pag = [float(r.vars['value2']) for r in result_set.members]
    Par = [float(r.vars['value3']) for r in result_set.members]
    Pdg = [float(r.vars['value4']) for r in result_set.members]
    Pdr = [float(r.vars['value5']) for r in result_set.members]
    gen = [float(r.vars['value6']) for r in result_set.members]
    PrepPtot = [] 
    for i in range(len(t)):
        Prep = Par[i] + Pdr[i]
        Ptot = Pag[i] + Par[i] + Pdg[i] + Pdr[i]
        PrepPtot.append(Prep/Ptot)
    axisC.plot(t, PrepPtot, 'k')  
    fs = 8 
axisC.set_xticks([0, 24, 48, 72])#, 96, 120, 144, 168, 192, 216, 240])
axisC.set_ylim([0, 0.45])
axisC.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
axisC.set_xlim([0, 72])
axisC.set_xlabel(r'Time (h)')
axisC.set_ylabel(r'$ \^\beta $ (colored), $P_{rep}$/$P_{tot}$ (black)')
axisC.set_title('Following Old Pole Cell Over Time')

pad = 0.04
side = "right"
divider = make_axes_locatable(axisC)
cax = divider.append_axes(side, size="25%", pad=pad, axisbg='none')
cax.set_xticklabels(['']*10)
cax.set_yticklabels(['']*10)
for spine in ['right', 'left', 'top', 'bottom']:
    cax.spines[spine].set_color('none')
cax.tick_params(top="off", bottom="off", right="off", left="off")

#Legend
x = [0.05, 0.2]
asnry = [0.9, 0.9]
cax.plot(x, asnry, color='red')
cax.text(0.25, 0.89, 'Asymmetric,', fontsize=fs)
cax.text(0.26, 0.845, 'No repair', fontsize=fs)
asory = [0.78, 0.78]
cax.plot(x, asory, color='#F781F3')
cax.text(0.25, 0.77, 'Asymmetric,', fontsize=fs)
cax.text(0.26, 0.725, ' Optimal repair', fontsize=fs)
asdry = [0.66, 0.66]
cax.plot(x, asdry, color='#fbbf07')
cax.text(0.25, 0.65, 'Asymmetric,', fontsize=fs)
cax.text(0.26, 0.61, ' Adaptive repair', fontsize=fs)
snry = [0.545, 0.545]
cax.plot(x, snry, color='blue')
cax.text(0.25, 0.535, 'Symmetric,', fontsize=fs)
cax.text(0.26, 0.49, 'No repair', fontsize=fs)
sory = [0.425, 0.425]
cax.plot(x, sory, color='#81BEF7')
cax.text(0.25, 0.415, 'Symmetric,', fontsize=fs)
cax.text(0.26, 0.37, 'Optimal repair', fontsize=fs)
sdry = [0.305, 0.305]
cax.plot(x, sdry, color='green')
cax.text(0.25, 0.295, 'Symmetric,', fontsize=fs)
cax.text(0.26, 0.25, 'Adaptive repair', fontsize=fs)
cax.set_xlim([0,1])
cax.set_ylim([0,1])



fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)