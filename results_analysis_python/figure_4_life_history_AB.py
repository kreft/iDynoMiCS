#!/usr/bin/python
from __future__ import division
import aging_extras
import os
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
input_path = os.path.join(base_path, 'results')
output_path = os.path.join(input_path, 'figure_4_life_history.pdf')

fig = toolbox_plotting.BmcFigure(double_column=True, height='double')

axisA = fig.add_subplot('A', 221)
paths = ['an10_life_history_results.xml', 'ao10_life_history_results.xml',
         'sn10_life_history_results.xml', 'so10_life_history_results.xml', 'T010AS12OBConstEnv0_80D_life_history_results.xml']
colorsA = ['red','#F781F3','blue','#81BEF7', '#fbbf07']
attributes1 = {'name':'specific growth rate', 'header':'time,value'}
attributes2 = {'name':'specific growth rate', 'name2':'age', 'header':'time,value,value2'}

xlim, ylim = 24, 0.6

for i in range(5):
    #attributes = attributes2    
    if i == 4:
        attributes = attributes2
    else:
        attributes = attributes1
    print attributes
    results_path = os.path.join(input_path, paths[i])
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        result.vars['time'] = float(result.vars['time'])
    result_set.members.sort(key=lambda r: r.vars['time'])
    t = [float(r.vars['time']) for r in result_set.members]
    #updated so as this takes value 1 - for my results this is specific growth
    #rate, which value 2 is age
    val = [float(r.vars['value']) for r in result_set.members]
    j = 0
    while j < len(t)-1:
        x = [t[j]]
        y = [val[j]]
        while (j < len(t)-1) and (abs(val[j+1]-val[j]) < 0.05):
            x.append(t[j+1])
            y.append(val[j+1])
            j += 1
        axisA.plot(x, y, color=colorsA[i])
        j += 1

fs = 8 

axisA.text(9, 0.36, 'Symmetric, No repair', color=colorsA[2], va='top',
                                                                   fontsize=fs)
axisA.text(9, 0.42, 'Symmetric, Optimal repair', color=colorsA[3], va='top',
                                                                   fontsize=fs)
axisA.text(1,0.56,'Asymmetric, No repair', color=colorsA[0], va='center',
                                                                   fontsize=fs)
axisA.text(2,0.48,'Asymmetric, Optimal repair', color=colorsA[1], va='center',
                                                                   fontsize=fs)
#these coordinates may need changing! This has removed the marker for 0 generation,
#so this will need replacing after this has been run
#axisA.text([1.5,0.51], 'Asymmetric, Dynamic repair', color=colorsA[4], va='center',
                                                                   #fontsize=fs)
#axisA.text(1.5,0.51,'0', fontsize=fs)
axisA.text(3.2, 0.395, '1', fontsize=fs)
axisA.text(3,0.25,'2',color=colorsA[0], fontsize=fs)
axisA.text(6,0.1,'3',color=colorsA[0], fontsize=fs)
axisA.text(4.5,0.32,'2',color=colorsA[1], fontsize=fs)
axisA.text(7,0.26,'3',color=colorsA[1], fontsize=fs)
axisA.text(10,0.21,'4',color=colorsA[1], fontsize=fs)
axisA.text(14,0.16,'5',color=colorsA[1], fontsize=fs)
axisA.text(18,0.11,'6',color=colorsA[1], fontsize=fs)
axisA.text(22,0.03,'7',color=colorsA[1], fontsize=fs)
aging_extras.draw_cell(axisA, 0.59*xlim, 0.69*xlim, 0.9*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisA, 0.79*xlim, 0.8*ylim, 0.99*xlim, ylim, y2xscale=ylim/xlim)
axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 4, 8, 12, 16, 20])#, 24])
axisA.set_ylim([0, ylim])
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
#This correlates to figure A and is to help the reader more than anything - these
    #numbers were by eye, not calculated
#axisB.text(0.021, 0.575,'0', color=colors[0], fontsize=fs)
#axisB.text(0.055, 0.36, '1', color=colors[0], fontsize=fs)
#axisB.text(0.019, 0.245, '2', color=colors[0], fontsize=fs)
#axisB.text(0.009, 0.08, '3', color=colors[0], fontsize=fs)
#axisB.text(0.009, 0.0, '4', color=colors[0], fontsize=fs)
#axisB.text(0.03, 0.515,'0', color=colors[1], fontsize=fs)
#axisB.text(0.075, 0.38, '1', color=colors[1], fontsize=fs)
#axisB.text(0.012, 0.35, '2', color=colors[1], fontsize=fs)
#axisB.text(0.0155, 0.31, '3', color=colors[1], fontsize=fs)
#axisB.text(0.006, 0.24, '4', color=colors[1], fontsize=fs)
#axisB.text(0.003, 0.16, '5', color=colors[1], fontsize=fs)
aging_extras.draw_cell(axisB, 0.6*xlim, 0.7*xlim, 0.9*ylim, 0.05*xlim,
                               y2xscale=ylim/xlim, toxic=True, arrow=0.01*xlim)
aging_extras.draw_const_env(axisB, 0.8*xlim, 0.8*ylim, xlim, ylim,
                                                            y2xscale=ylim/xlim)
axisB.set_xlim([0.00, xlim])
axisB.set_xticks([0.00, 0.02, 0.04, 0.06, 0.08])
axisB.set_ylim([0, ylim])
axisB.set_xlabel('Frequency in population')
axisB.set_title('Growth Rate Distribution')
#Only needed when no A axis - also make y tick labels true if A absent
#axisB.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
setp( axisB.get_yticklabels(), visible=False)

fig.process_subplots()
#axisB.tick_params(left='off')
#axisD.tick_params(left='off')
#axisCD.tick_params(top="off", bottom="off", right="off", left="off")
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)


