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


axisA.set_xlim([0, xlim])
axisA.set_xticks([0, 24, 48, 72])
axisA.set_ylim([0, ylim])
axisA.set_xlabel(r'Time (h)')
axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)')
axisA.set_title('Following Old Pole Cell Over Time')


fig.process_subplots()
#axisB.tick_params(left='off')
#axisD.tick_params(left='off')
#axisCD.tick_params(top="off", bottom="off", right="off", left="off")
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)


