from __future__ import division
import os
from pylab import *
import toolbox_basic as basic
import toolbox_results as results
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy
import math
import aging_extras
import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable
from statsmodels.stats.proportion import proportions_ztest
fs, fs2, fs3 = 10, 10, 8
stype = '.png'

constant = ['constant_SAR_ASNR(20170504_2048)', 'constant_SAR_ASNR(20170723_0023)', 'constant_SAR_ASNR(20170514_1151)', 'constant_SAR_ASNR(20170716_1621)', 'constant_SAR_ASNR(20170716_2311)', 'constant_SAR_ASNR(20170907_1731)', 'constant_SAR_ASNR(20170907_2341)', 'constant_SAR_ASNR(20170912_1214)', 'constant_SAR_ASNR(20171011_0955)', 'constant_SAR_ASNR(20171017_1221)',
            'constant_SOB_ASNR(20170505_1238)', 'constant_SOB_ASNR(20170510_2246)', 'constant_SOB_ASNR(20170514_1703)', 'constant_SOB_ASNR(20170716_1816)', 'constant_SOB_ASNR(20170717_1153)', 'constant_SOB_ASNR(20170907_1944)', 'constant_SOB_ASNR(20170908_0132)', 'constant_SOB_ASNR(20170913_0747)', 'constant_SOB_ASNR(20171011_1722)', 'constant_SOB_ASNR(20171011_1737)',
            'constant_SAR_SOB(20170505_0912)', 'constant_SAR_SOB(20170723_0124)', 'constant_SAR_SOB(20170514_1219)', 'constant_SAR_SOB(20170716_1657)', 'constant_SAR_SOB(20170716_2340)', 'constant_SAR_SOB(20170907_1756)', 'constant_SAR_SOB(20170908_0010)', 'constant_SAR_SOB(20170912_1432)', 'constant_SAR_SOB(20171011_1019)', 'constant_SAR_SOB(20171011_1020)']
dynamic = [ 'dynamic_SAR_ASNR(20170505_1305)', 'dynamic_SAR_ASNR(20170723_0242)', 'dynamic_SAR_ASNR(20170514_1151)', 'dynamic_SAR_ASNR(20170716_2143)', 'dynamic_SAR_ASNR(20170717_1226)', 'dynamic_SAR_ASNR(20170907_2005)', 'dynamic_SAR_ASNR(20170908_0152)', 'dynamic_SAR_ASNR(20170913_1013)', 'dynamic_SAR_ASNR(20171011_1809)', 'dynamic_SAR_ASNR(20171011_1905)',
           'dynamic_SOB_ASNR(20170510_2233)', 'dynamic_SOB_ASNR(20170511_1456)', 'dynamic_SOB_ASNR(20170515_0938)', 'dynamic_SOB_ASNR(20170719_1921)', 'dynamic_SOB_ASNR(20170723_0326)', 'dynamic_SOB_ASNR(20170911_1620)', 'dynamic_SOB_ASNR(20170916_0207)', 'dynamic_SOB_ASNR(20170923_0723)', 'dynamic_SOB_ASNR(20171020_2334)', 'dynamic_SOB_ASNR(20171021_0454)',
            'dynamic_SAR_SOB(20170505_1355)', 'dynamic_SAR_SOB(20170510_2247)', 'dynamic_SAR_SOB(20170514_1239)', 'dynamic_SAR_SOB(20170716_2235)', 'dynamic_SAR_SOB(20170717_1343)', 'dynamic_SAR_SOB(20170907_2055)', 'dynamic_SAR_SOB(20170908_0240)', 'dynamic_SAR_SOB(20170913_1517)', 'dynamic_SAR_SOB(20171017_1221)', 'dynamic_SAR_SOB(20171017_1321)']

title1, title2, title3 = ['AR/DS', 'FR/DS', 'AR/FR']
titles = ['AR/DS', 'FR/DS', 'AR/FR']

def get_biomass(sim_dir_path, a):
    times1, biomass1, population1, growthrate1 = [], [], [], []
    times2, biomass2, population2, growthrate2 = [], [], [], []
    file_path = '/Volumes/Robyn_W_2/july_2018/paper/'+sim_dir_path
    file_dir = file_path+'/agent_sum'
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        species1 = results.SpeciesOutput(output, 'OldieA')
        requirements={}
        cells = species1.find_cells(requirements)
        times1.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population1.append(float(cell.vars['population']))
        biomass1.append(float(cell.vars['mass']))
        growthrate1.append(float(cell.vars['growthRate']))
        species2 = results.SpeciesOutput(output, 'OldieB')
        cells = species2.find_cells(requirements)
        times2.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population2.append(float(cell.vars['population']))
        biomass2.append(float(cell.vars['mass']))
        growthrate2.append(float(cell.vars['growthRate']))
    basic.rm_dir(file_dir)
    lists_a, lists_b = [population1, biomass1, growthrate1], [population2, biomass2, growthrate2]
    t_a1, t_b1 = times1, times2
    for i in range(3):
        list1, list2, list3, list4 = t_a1, lists_a[i], t_b1, lists_b[i]
        t_a, lists_a[i] = (list(x) for x in zip(*sorted(zip(list1, list2), key=itemgetter(0))))
        t_b, lists_b[i] = (list(x) for x in zip(*sorted(zip(list3, list4), key=itemgetter(0))))
    biomass_ratio, population_ratio, growthrate_ratio = [], [], []
    for j in range(len(lists_a[1])):
        if lists_a[1][j] == 0 and lists_b[1][j] == 0 or lists_a[1][j] == 0 and lists_b[1][j] != 0 or lists_a[1][j] != 0 and lists_b[1][j] == 0:
            biomass_ratio.append(0)
        else:
            biomass_ratio.append(numpy.log(lists_a[1][j]/lists_b[1][j]))
        if lists_a[0][j] == 0 and lists_b[0][j] == 0 or lists_a[0][j] == 0 and lists_b[0][j] != 0 or lists_a[0][j] != 0 and lists_b[0][j] == 0:
            population_ratio.append(0)
        else:
            population_ratio.append(lists_a[0][j]/lists_b[0][j])
        if lists_a[2][j] == 0 and lists_b[2][j] == 0  or lists_a[2][j] == 0 and lists_b[2][j] != 0 or lists_a[2][j] != 0 and lists_b[2][j] == 0:
            growthrate_ratio.append(0)
        else:
            growthrate_ratio.append(lists_a[2][j]/lists_b[2][j])
    return t_a, biomass_ratio, population_ratio, growthrate_ratio

counts = [[0, 0, 0], [0, 0, 0]]
def plot_all(condyn, plotting=0, axis=0, ylabel='biomass ratio'):
    string = '23'
    string1 = axis+1
    string2, string3 = string1+1, string1+2
    
    caxB = fig.add_subplot(334)
    caxB.set_xticklabels(['']*10)
    caxB.set_yticklabels(['']*10)
    caxB.spines['left'].set_color('none')
    caxB.spines['right'].set_color('none')
    caxB.spines['top'].set_color('none')
    caxB.spines['bottom'].set_color('none')
    caxB.tick_params(top="off", right="off", bottom="off", left="off")
    caxB.set_ylabel(ylabel+'\n \n \n')
    axisA = fig.add_subplot(int(string+str(string1)))
    axisB, axisC = fig.add_subplot(int(string+str(string2)), sharey=axisA), fig.add_subplot(int(string+str(string3)), sharey=axisA)
    axis = [axisA, axisB, axisC]
    if condyn == constant:
        axisA.text(.05, .9, 'A', horizontalalignment='left', transform=axisA.transAxes)
        axisB.text(.05, .9, 'B', horizontalalignment='left', transform=axisB.transAxes)
        axisC.text(.05, .9, 'C', horizontalalignment='left', transform=axisC.transAxes)
        axisA.set_title(title1)
        axisB.set_title(title2)
        axisC.set_title(title3)
        axisA.set_ylabel('Constant environment', fontsize=10)
        axisC.set_xlim([0, 2100])
        prefix = 'constant_competitions/'
    elif condyn == dynamic:
        axisA.text(.05, .9, 'D', horizontalalignment='left', transform=axisA.transAxes)
        axisB.text(.05, .9, 'E', horizontalalignment='left', transform=axisB.transAxes)
        axisC.text(.05, .9, 'F', horizontalalignment='left', transform=axisC.transAxes)
        axisA.set_ylabel('Chemostat environment', fontsize=10)
        axisA.set_xlabel('Time (h)', fontsize=10)
        axisB.set_xlabel('Time (h)', fontsize=10)
        axisC.set_xlabel('Time (h)', fontsize=10)
        prefix = 'dynamic_competitions/'
    setp(axisB.get_yticklabels(), visible=False)
    setp(axisC.get_yticklabels(), visible=False)
    if condyn==dynamic:
        axisC.set_xticks([0, 2000, 4000, 6000, 8000])
        axisA.set_xticks([0, 50, 100, 150, 200])
        axisC.set_xlim([0, 8500])
        coun = 1
    else:
        coun = 0
    for a in range(30):
        if a < 10: 
            axiss = axis[0]
            app = 0
        elif a < 20: 
            axiss = axis[1]
            app = 1
        else: 
            axiss = axis[2]
            app = 2
        time, biomass, population, growthrate = get_biomass(prefix+condyn[a], a)
        ratios = [biomass, population, growthrate]
        px, py = [], []
        for b in range(len(time)):
            if b == len(time)-1:
                continue
            px.append(time[b])
            py.append(ratios[plotting][b])
        axiss.plot(px, py, linewidth=1.5)
        if py[-1] > 0:
            counts[coun][app] += 1
        axiss.plot([time[0], time[len(time)-1]], [0, 0], 'k--')
    return

fig = plt.figure(figsize=(8.27, 5.51))
plot_all(constant, 0, 0, 'Log biomass ratio')
plot_all(dynamic, 0, 3, 'Log biomass ratio')

fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.2, hspace=0.2)
os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/const_chemo_comps/')
plt.savefig('All ratios'+stype, bbox_inches='tight', dpi=600)
plt.close()

printing = ['Constant', 'Dynamic']
stats = []
for a in range(len(counts)):
    for b in range(len(counts[a])):
        stat, pval = proportions_ztest(int(counts[a][b]), 10, 0.5)
        stats.append([printing[a], titles[b], stat, pval])
with open('Constant and chemostat competition stats.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Environment', 'Competition', 'Statistic', 'P value'])
    for a in stats:
        writer.writerow(a)
