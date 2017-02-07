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
input_path = os.path.join(base_path, 'results', 'comparison_clegg', 'Final_versions', 'Reproducing_figure_4')
output_path = os.path.join(input_path, 'figure_4_life_history_ABC_3D_adapt.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')

axisA = fig.add_subplot('A', 221)
paths = ['T010S12OBConstEnvCleggComp_life_history_results.xml', 'T010S12NRConstEnvCleggComp_life_history_results.xml', 'T010S12007ConstEnvCleggComp_life_history_results.xml',
         'T010AS12NRConstEnvCleggComp_life_history_results.xml', 'T010AS12007ConstEnvCleggComp_life_history_results.xml', 'T010AS12OBConstEnv1_0D_Ind_life_history_results.xml']
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
#these coordinates may need changing! This has removed the marker for 0 generation,
#so this will need replacing after this has been run
#axisA.text([1.5,0.51], 'Asymmetric, Dynamic repair', color=colorsA[4], va='center',
                                                                   #fontsize=fs)
axisA.text(1.25,0.55,'0', color='#8B4513', fontsize=fs)
axisA.text(0.8, 0.425, '1', color='#8B4513', fontsize=fs)
axisA.text(3.15,0.3,'2', color='#8B4513', fontsize=fs)
axisA.text(5,0.15,'3',color=colorsA[3], fontsize=fs)
axisA.text(8,0,'4',color=colorsA[3], fontsize=fs)
axisA.text(5.75,0.225,'3', color='#8B4513', fontsize=fs)
axisA.text(9,0.225,'4', color='#8B4513', fontsize=fs)
axisA.text(12.5,0.125,'5',color=colorsA[4], fontsize=fs)
axisA.text(19.5,0.08,'6',color=colorsA[4], fontsize=fs)
axisA.text(22,0.04,'7',color=colorsA[4], fontsize=fs)
axisA.text(12,0.235,'5',color=colorsA[5], fontsize=fs)
axisA.text(15,0.235,'6',color=colorsA[5], fontsize=fs)
axisA.text(17.5,0.235,'7',color=colorsA[5], fontsize=fs)
axisA.text(20,0.235,'8',color=colorsA[5], fontsize=fs)
axisA.text(23,0.235,'9',color=colorsA[5], fontsize=fs)
aging_extras.draw_cell(axisA, 0.59*xlim, 0.69*xlim, 0.9*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisA, 0.79*xlim, 0.8*ylim, 0.99*xlim, ylim, y2xscale=ylim/xlim)
axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 4, 8, 12, 16, 20, 24])
#axisA.set_ylim([0, ylim])
axisA.set_xlabel(r'Time (h)')
axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
axisA.set_title('Following Old Pole Cell Over Time')



axisB = fig.add_subplot('B', 222)
paths = ['T010AS12NR100_results.xml', 'T010AS12007100_results.xml', 'T010AS12OB100_results.xml']
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

axisB.text(0.026,0.58,'0', color=colorsB[2], fontsize=fs)
axisB.text(0.041,0.51,'1', color=colorsB[2], fontsize=fs)
axisB.text(0.065,0.38,'2', color=colorsB[2], fontsize=fs)
axisB.text(0.025,0.27,'3-', color=colorsB[2], fontsize=fs)
axisB.text(0.022,0.58,'0', color=colorsB[0], fontsize=fs)
axisB.text(0.038,0.51,'1', color=colorsB[0], fontsize=fs)
axisB.text(0.061,0.38,'2', color=colorsB[0], fontsize=fs)
axisB.text(0.025,0.24,'3', color=colorsB[0], fontsize=fs)
axisB.text(0.014,0.08,'4', color=colorsB[0], fontsize=fs)
axisB.text(0.008,0.00,'5', color=colorsB[0], fontsize=fs)
axisB.text(0.033,0.54,'0', color=colorsB[1], fontsize=fs)
axisB.text(0.045,0.48,'1', color=colorsB[1], fontsize=fs)
axisB.text(0.090,0.40,'2', color=colorsB[1], fontsize=fs)
axisB.text(0.021,0.33,'3', color=colorsB[1], fontsize=fs)

aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
axisB.set_xlim([0, 0.105])
axisB.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10])
axisB.set_ylim([0, 0.6])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate: Distribution')
#Only needed when no A axis - also make y tick labels true if A absent
#axisB.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
setp( axisB.get_yticklabels(), visible=False)

pad = 0.00
side = "right"
divider = make_axes_locatable(axisB)
caxB = divider.append_axes(side, size="25%", pad=pad, axisbg='none', frameon=True)
caxB.set_xticklabels(['']*10)
caxB.set_yticklabels(['']*10)
caxB.spines['left'].set_color('none')
caxB.tick_params(top="off", right="off", bottom="off", left="off")
caxB.scatter(0.5, 0.403227, color=colorsB[0])
plt = ([0.5,0.5])
NR = ([0.403227+0.148353, 0.403227-0.148353])
caxB.plot(plt, NR, '--', color=colorsB[0])
caxB.plot([0.45, 0.55], [NR[0], NR[0]], color=colorsB[0])
caxB.plot([0.45, 0.55], [NR[1], NR[1]], color=colorsB[0])
caxB.scatter(1, 0.422097, color=colorsB[1])
plt = ([1,1])
opt = ([0.422097+0.095199, 0.422097-0.095199])
caxB.plot(plt, opt, '--', color=colorsB[1])
caxB.plot([0.95, 1.05], [opt[0], opt[0]], color=colorsB[1])
caxB.plot([0.95, 1.05], [opt[1], opt[1]], color=colorsB[1])
caxB.scatter(1.5, 0.424253, color=colorsB[2])
plt = ([1.5, 1.5])
AR = ([0.424253+0.110226, 0.424253-0.110226])
caxB.plot(plt, AR, '--', color=colorsB[2])
caxB.plot([1.45, 1.55], [AR[0], AR[0]], color=colorsB[2])
caxB.plot([1.45, 1.55], [AR[1], AR[1]], color=colorsB[2])
caxB.set_xlim([0, 2])
caxB.set_ylim([0, 0.6])
caxB.set_title(r'Means $\pm$ SD')
x = ([0.5, 1, 1.5])
my_xticks = (['No repair', 'Optimal repair', 'Adaptive repair'])


axisC = fig.add_subplot('C', 212)
paths = ['T010AS12OBConstEnv10D_Ind_life_history_results.xml', 'T010S12OB3day_life_history_results.xml']
pathsprot = ['T010AS12OBConstEnv10D_Ind_results_for_Prep_Ptot.xml', 'T010S12OB3day_results_for_Prep_Ptot.xml']
colorsC = ['#fbbf07', 'green', '#81BEF7']
attributes = {'name':'specific growth rate', 'name2':'age', 'name3':'generation', 'header':'time,value,value2,value3'}
attributesprot = {'name':'age', 'name2':'activeBiomassGrowth', 'name3':'activeBiomassRepair', 'name4':'inactiveBiomassGrowth',
'name5':'inactiveBiomassRepair', 'name6':'generation', 'header':'time,value,value2,value3,value4,value5,value6'}

muS = 0.6
Yr = 0.8

for m in range(2):
    results_path = os.path.join(input_path, pathsprot[m])
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
    for n in range(len(t)):
        Prep = Par[n] + Pdr[n]
        Ptot = Pag[n] + Par[n] + Pdg[n] + Pdr[n]
        PrepPtot.append(Prep/Ptot)
    axisC.plot(t, PrepPtot, 'k')  
    fs = 8 

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
    j = 0
    while j < len(t)-1:
        x, y = [], []
        while (j < len(t)-1) and [gen[j]] == [gen[j+1]]:
            if t[j] < 72:
                age = val[j+1]
                x.append(t[j+1])
                y.append([float(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))])
            j += 1
        axisC.plot(x, y, color=colorsC[i])
        j += 1       
    fs = 8
xall = [0, 72]
yall = [0.07, 0.07]
axisC.plot(xall, yall, color=colorsC[2])
    

axisC.set_xticks([0, 24, 48, 72])#, 96, 120, 144, 168, 192, 216, 240])
axisC.set_ylim([0, 0.45])
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