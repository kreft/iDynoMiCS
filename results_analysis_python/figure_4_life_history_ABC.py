#!/usr/bin/python
from __future__ import division
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_plotting
import toolbox_results

#base_path = os.path.join('~', 'Dropbox', 'EclipseWorkspace', 'iDynoAge')
#base_path = os.path.join('C:\\', 'Users', 'RJW598', 'git', 'iDynoMiCS')
base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'Figure_4')
output_path = os.path.join(input_path, 'figure_4_life_history_ABC.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')

axisA = fig.add_subplot('A', 221)
paths = ['T010AS12NRConstEnvCleggComp_life_history_results.xml', 'T010AS12007ConstEnvCleggComp_life_history_results.xml',
         'T010S12NRConstEnvCleggComp_life_history_results.xml', 'T010S12007ConstEnvCleggComp_life_history_results.xml', 'T010AS12OBConstEnv1_0D_Ind_life_history_results.xml', 'T010S12OBConstEnvCleggComp_life_history_results.xml']
colorsA = ['red','#F781F3','blue','#81BEF7', '#fbbf07', 'green']
#attributes1 = {'name':'specific growth rate', 'header':'time,value'}
attributes = {'name':'specific growth rate', 'name2':'age', 'name3':'generation', 'header':'time,value,value2,value3'}

xlim, ylim = 25, 0.6

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
axisA.text(5,0.15,'3',color=colorsA[0], fontsize=fs)
axisA.text(8,0,'4',color=colorsA[0], fontsize=fs)
axisA.text(5.75,0.225,'3', color='#8B4513', fontsize=fs)
axisA.text(9,0.225,'4', color='#8B4513', fontsize=fs)
axisA.text(12.5,0.125,'5',color=colorsA[1], fontsize=fs)
axisA.text(19.5,0.08,'6',color=colorsA[1], fontsize=fs)
axisA.text(22,0.04,'7',color=colorsA[1], fontsize=fs)
axisA.text(12,0.235,'5',color=colorsA[4], fontsize=fs)
axisA.text(15,0.235,'6',color=colorsA[4], fontsize=fs)
axisA.text(17.5,0.235,'7',color=colorsA[4], fontsize=fs)
axisA.text(20,0.235,'8',color=colorsA[4], fontsize=fs)
axisA.text(23,0.235,'9',color=colorsA[4], fontsize=fs)
aging_extras.draw_cell(axisA, 0.59*xlim, 0.69*xlim, 0.9*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisA, 0.79*xlim, 0.8*ylim, 0.99*xlim, ylim, y2xscale=ylim/xlim)
axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 4, 8, 12, 16, 20, 24])
#axisA.set_ylim([0, ylim])
axisA.set_xlabel(r'Time (h)')
axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
axisA.set_title('Following Old Pole Cell Over Time')



axisB = fig.add_subplot('B', 222)
paths = ['an10_const_steady_results.xml', 'ao10_const_steady_results.xml', 'T010AS12OBConstEnvH110D_results.xml']
attributes = {'name':'specific growth rate population structure',
                    'starting_time':'2400', 'header':'bin,frequency'}
colorsB = ['red','#F781F3','#fbbf07']

xlim, ylim = 0.082, 0.6

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

axisB.text(0.025,0.55,'0', color=colorsA[0], fontsize=fs)
axisB.text(0.03,0.52,'0', color=colorsA[4], fontsize=fs)
axisB.text(0.04,0.48,'0', color=colorsA[1], fontsize=fs)
axisB.text(0.06,0.38,'1', color=colorsA[0], fontsize=fs)
axisB.text(0.06,0.43,'1', color=colorsA[4], fontsize=fs)
axisB.text(0.082,0.4,'1', color=colorsA[1], fontsize=fs)
aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
axisB.set_xlim([0, 0.09])
axisB.set_xticks([0, 0.02, 0.04, 0.06, 0.08])
axisB.set_ylim([0, ylim])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate Distribution')
#Only needed when no A axis - also make y tick labels true if A absent
#axisB.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
setp( axisB.get_yticklabels(), visible=False)


axisC = fig.add_subplot('C', 223)
paths = ['T010AS12OBConstEnv1_0D_Ind_life_history_results.xml', 'T010S12OBConstEnvCleggComp_life_history_results.xml']
colorsC = ['#fbbf07', 'green']
attributes = {'name':'specific growth rate', 'name2':'age', 'name3':'generation', 'header':'time,value,value2,value3'}

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
    j = 0
    while j < len(t)-1:
        x, y = [], []
        while (j < len(t)-1) and [gen[j]] == [gen[j+1]]:    
            age = val[j+1]
            x.append(t[j+1])
            y.append([float(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))])
            j += 1
        axisC.plot(x, y, color=colorsC[i])
        j += 1       
    fs = 8 
axisC.set_xticks([0, 4, 8, 12, 16, 20, 24])
axisC.text(1,0.05,'0',color=colorsA[4], fontsize=fs)
axisC.text(3.75,0.1,'1',color=colorsA[4], fontsize=fs)
axisC.text(6,0.2,'2',color=colorsA[4], fontsize=fs)
axisC.text(9,0.2,'3',color=colorsA[4], fontsize=fs)
axisC.text(11.5,0.18,'4',color=colorsA[4], fontsize=fs)
axisC.text(15,0.15,'5',color=colorsA[4], fontsize=fs)
axisC.text(17.5,0.15,'6',color=colorsA[4], fontsize=fs)
axisC.text(20.25,0.13,'7',color=colorsA[4], fontsize=fs)
axisC.text(22,0.2,'8',color=colorsA[4], fontsize=fs)
axisC.text(24.5,0.27,'9',color=colorsA[4], fontsize=fs)
axisC.set_xlabel(r'Time (h)')
axisC.set_ylabel(r'Dynamic beta')
axisC.set_title('Following Old Pole Cell Over Time')

#move this to the bottom if you actually want there to be a sub plot D!
fig.process_subplots()
#Legend
axisD = fig.add_subplot('D', 224, frameon=False)
x = [0.1, 0.2]
asnry = [0.9, 0.9]
axisD.plot(x, asnry, color='red')
axisD.text(0.25, 0.9, 'Asymmetric, No repair', fontsize=fs)
asory = [0.825, 0.825]
axisD.plot(x, asory, color='#F781F3')
axisD.text(0.25, 0.825, 'Asymmetric, Optimal repair', fontsize=fs)
asdry = [0.75, 0.75]
axisD.plot(x, asdry, color='#fbbf07')
axisD.text(0.25, 0.75, 'Asymmetric, Dynamic repair', fontsize=fs)
snry = [0.675, 0.675]
axisD.plot(x, snry, color='blue')
axisD.text(0.25, 0.675, 'Symmetric, No repair', fontsize=fs)
sory = [0.6, 0.6]
axisD.plot(x, sory, color='#81BEF7')
axisD.text(0.25, 0.6, 'Symmetric, Optimal repair', fontsize=fs)
sdry = [0.525, 0.525]
axisD.plot(x, sdry, color='green')
axisD.text(0.25, 0.525, 'Symmetric, Dynamic repair', fontsize=fs)

axisD.set_xlim([0.00, 1.00])
axisD.set_ylim([0.00, 1.00])


fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)