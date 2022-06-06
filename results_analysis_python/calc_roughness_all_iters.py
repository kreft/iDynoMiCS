#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import aging_extras
import os
import math
from pylab import *
import toolbox_basic
import toolbox_results
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting_age as toolbox_plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import numpy
import matplotlib.pyplot as plt
import calc_roughness

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = toolbox_basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'roughness')
output_path = os.path.join(input_path, 'roughness.pdf')

#fig = toolbox_plotting.BmcFigure(double_column=True) #, height='double')
colors = ['r','b','g']

#axisABC = fig.add_subplot('', 111, frameon=True)
high = ['NANRH154_25dayhalf(20160301_1158)', 'NANRH152_25day(20160315_0935)', 'NANRH152_25day(20160317_1147)', 'NANRH152_25day(20160321_1451)']
medium = ['NANRM154_25dayhalf(20160301_2312)', 'NANRM152_25day(20160316_1440)', 'NANRM152_25day(20160318_1835)', 'NANRM152_25day(20160322_1524)']
low = ['NANRL154_25dayhalf(20160301_1343)', 'NANRL152_25day(20160315_1337)', 'NANRL152_25day(20160317_1534)', 'NANRL152_25day(20160321_1746)']
for i in range(1):
    if i==0: paths = low
    if i==1: paths = medium
    if i==3: paths = low
    heightgrid1 = numpy.zeros((256/4, 256), dtype=float)
    heightgrid2 = numpy.zeros((256/4, 256), dtype=float)
    heightgrid3 = numpy.zeros((256/4, 256), dtype=float)
    heightgrid4 = numpy.zeros((256/4, 256), dtype=float)
    heightgridall = numpy.zeros((256/4, 256*4), dtype=float)
    heightmean = []
    heightstd = []
    heighttoplot = []
    biofilmheight = []
    heightstdtoplot = []
    for j in range(4):
        roughnessj, heightj = roughness[j], height[j]
        sim1 = os.path.join(input_path, paths[j])
        sim = toolbox_idynomics.SimulationDirectory(sim1)
        print sim
        iter_info = sim.get_last_iterate_number()        
        alliters = sim.get_iterate_numbers()
        if j==0: heightgrid=heightgrid1
        if j==1: heightgrid=heightgrid2
        if j==2: heightgrid=heightgrid2
        if j==3: heightgrid=heightgrid4
        for k in range(len(alliters)):
            calc_roughness.calc_roughness(sim1, alliters[k])
            '''            
            high = int(high)
            if j==0: heightgrid1[high, k]=rough
            if j==1: heightgrid2[high,k]=rough
            if j==2: heightgrid3[high,k]=rough
            if j==3: heightgrid4[high,k]=rough
            '''
        sim.clean_up()
    '''
    for a in range(65):
        for c in range(256):
            heightgridall[a,c]=heightgrid1[a,c]
            heightgridall[a,(c+256)]=heightgrid2[a,c]
            heightgridall[a,(c+512)]=heightgrid2[a,c]
            heightgridall[a,(c+768)]=heightgrid2[a,c]
    binsofrows = []
    for d in range(65):
        row = heightgridall[d,:]
        newrow = []
        for e in range(len(row)):
            if not row[e]==0:
                if row[e]=='NaN':
                    row[e] = row[e-1]
                newrow.append(row[e])
        mean = numpy.mean(newrow)
        std = numpy.std(newrow)
        heightmean.append(mean)
        heightstd.append(std)
    allheights = range(0, len(heightmean))
    axisABC.plot(heightmean, allheights)
axisABC.set_xticks([0, 1, 2, 3, 4, 5, 6, 7])
axisABC.set_ylim([0, 256])
fig.process_subplots()
fig.subplots_adjust(left=0.09,right=0.99,top=0.96,bottom=0.08,wspace=0.0, hspace=0.35)
fig.save(output_path)
'''