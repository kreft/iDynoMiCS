from __future__ import division
import os
#from pylab import *
import toolbox_basic as basic
import toolbox_results as results
import toolbox_results
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from scipy.stats.mstats import gmean
#import operator
import aging_extras
#from geoMean import geomMeanExtension
#import geostdmod as geostd
#import struct

stype = '.png'

def build_life_history(sim_dir_path, attribute1='specific growth rate', attribute2='age', cell_id=(1,0),
                        biomass_names=['activeBiomassGrowth', 'activeBiomassRepair',
                                       'inactiveBiomassGrowth', 'inactiveBiomassRepair']):
    
    sim_dir_path = os.path.join('//Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/'+sim_dir_path+'/')
    time, generation, spec_growth, age, inact_rep, act_rep, inact_gro, act_gro, total_biomass = [], [], [], [], [], [], [], [], []
    requirements = {'family':str(cell_id[0]), 'genealogy':str(cell_id[1])}
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        #print filename
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output)
        cells = species.find_cells(requirements)
        time.append(output.time)
        if len(cells) == 0:
            continue
        cell = cells[0]
        generation.append(cell.vars['generation'])
        ar = cell.vars['activeBiomassRepair']
        act_rep.append(ar)
        ir = cell.vars['inactiveBiomassRepair']
        inact_rep.append(ir)
        ag = cell.vars['activeBiomassGrowth']
        act_gro.append(ag)
        ig = cell.vars['inactiveBiomassGrowth']
        inact_gro.append(ig)
        total_biomass.append(ar+ir+ag+ig)
        spec_growth.append(cell.get_specific_growth_rate(biomass_names))
        age.append(cell.vars[attribute2])
    basic.rm_dir(file_dir)
    return time, generation, spec_growth, age, inact_rep, act_rep, inact_gro, act_gro

def get_population(sim_dir_path):
    sim_dir_path = os.path.join('/Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/'+sim_dir_path+'/lastIter/')
    output = results.AgentOutput(path=sim_dir_path+'agent_State(last).xml')
    species = results.SpeciesOutput(output)
    cells = species.findcells()
    age, total_biomass = [], []
    for cell in cells:
        if len(cells) == 0:
            continue
        ar = cell.vars['activeBiomassRepair']
        ir = cell.vars['inactiveBiomassRepair']
        ag = cell.vars['activeBiomassGrowth']
        ig = cell.vars['inactiveBiomassGrowth']
        total_biomass.append(ar+ir+ag+ig)
        age.append(cell.vars['age'])
    return age, biomass

def build_population_structure(sim_dir_path, attribute='specific growth rate', 
                        bins=numpy.linspace(0, 0.6, 90), starting_time=2400,
                        biomass_names=['activeBiomassGrowth', 'activeBiomassRepair',
                                       'inactiveBiomassGrowth', 'inactiveBiomassRepair']):
    sim_dir_path = os.path.join('/Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/'+sim_dir_path+'/')
    attr_values = []
    total_pop = 0.0
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    all_growth = []
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        if output.time >= starting_time:
            species = results.SpeciesOutput(output)
            species.set_biomass_names(biomass_names)
            attr_values.extend(species.get_attribute_values(attribute))
            all_growth.append(species.get_attribute_values(attribute))
            total_pop += species.population()
    hist, bin_edges = numpy.histogram(attr_values, bins)
    bin_using, frequencies = [], []
    for i in range(len(hist)):
        bin_using.append((bins[i+1]+bins[i])/2)
        frequencies.append(float(hist[i]/total_pop))
    basic.rm_dir(file_dir)
    return bin_using, frequencies, all_growth

names2 = []
names = sorted(os.listdir('/Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/old_pole/'))
names3 = sorted(os.listdir('/Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/steady_state/'))
for a in names3:
    if a[13:15] == 'AS':
        names2.append(a)


ages, biomasses = [], []
colors = ['red', 'blue', '#F781F3', '#81BEF7', '#fbbf07', 'green']
colorsA = [colors[0], colors[2]]
colorsB = []
colorsC = []
ms, fs = 1.25, 6
for a in range(3):
    fig = plt.figure(figsize=(2.5, 8.27))
    axisA, axisB, axisC = fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)
    axisA.set_title('A', loc='left'), axisB.set_title('B', loc='left'), axisC.set_title('C', loc='left')
    ax = [axisC, axisA, axisB, axisC, axisA, axisB]
    color = [colors[4], colors[0], colors[2], colors[5], colors[1], colors[3]]
    #axisA, axisA, axisB, axisB, axisC, axisC]
    #color = [colors[0], colors[1], colors[2], colors[3], colors[4], colors[5]]
    for b in range(6):
        c = (b*3)+a
        input_path = os.path.join('/Volumes/Robyn_W_2/july_2018/paper/comparison_clegg/steady_state/'+names3[c]+'/lastIter/')
        results_path = os.path.join(input_path, 'agent_State(last).xml')
        agent_output = toolbox_results.AgentOutput(path=results_path)
        species = toolbox_results.SpeciesOutput(agent_output)
        aging_extras.scatter_population(ax[b], species, 'totalBiomass', 'age', color=color[b], markersize=ms)
    axisA.set_xlim([285, 630])
    axisA.set_ylim([0.0, 1.0])
    axisA.set_title('Age & Size Distributions', fontsize=10)
    plt.setp(axisA.get_xticklabels(), visible=False)
    axisB.set_xlim([285, 630])
    axisB.set_ylim([0.0, 1.0])
    axisB.set_ylabel(r'Cellular age $P_{dam}/(P_{act}+P_{dam})$', fontsize=10)
    plt.setp(axisB.get_xticklabels(), visible=False)
    axisC.set_xlim([285, 630])
    axisC.set_ylim([0.0, 1.0])
    axisC.set_xticks([300, 400, 500, 600])
    axisC.set_xlabel('Cellular total biomass (fg)', fontsize=10)
    left, right = 500, 620
    top, bottom = 0.98, 0.78
    asym_y = 0.915
    sym_y = 0.825
    axisA.text(290, 0.10, 'Generation 0', color=colors[0], fontsize=fs, ha='left', va='top')
    axisA.text(290, 0.3, '1', color=colors[0], fontsize=fs, ha='left', va='top')
    axisA.text(290, 0.65, '2', color=colors[0], fontsize=fs, ha='left', va='top')
    axisA.text(290, 0.86, '3', color=colors[0], fontsize=fs, ha='left', va='top')
    #axisA.text(290, 0.97, '4', color=colors[0], fontsize=fs, ha='left', va='top')
    axisA.plot([left+20, right, right, left+20, left+20], [bottom, bottom, top, top, bottom], '0.5')
    axisA.plot([left + 10 +20], [asym_y], 'o', color=colors[0], markeredgecolor='none', markersize=3)
    axisA.text(left + 20 +20, asym_y, 'Asymmetric, \n no repair', va='center', ha='left', fontsize=fs)
    axisA.plot([left + 10 +20], [sym_y], 'o', color=colors[1], markeredgecolor='none', markersize=3)
    axisA.text(left + 20 +20, sym_y, 'Symmetric, \n no repair', va='center', ha='left', fontsize=fs)
    axisA.set_title('Age & Size Distributions', fontsize=9)
    axisB.text(290, 0.1, 'Generation 0', color=colors[2], fontsize=fs, ha='left', va='top')
    axisB.text(290, 0.27, '1', color=colors[2], fontsize=fs, ha='left', va='top')
    axisB.text(290, 0.45, '2', color=colors[2], fontsize=fs, ha='left', va='top')
    axisB.text(290, 0.65, '3', color=colors[2], fontsize=fs, ha='left', va='top')
    axisB.text(290, 0.90, '4', color=colors[2], fontsize=fs, ha='left', va='top')
    #axisB.text(290, 0.75, '5', color=colors[2], fontsize=fs, ha='left', va='top')
    #axisB.text(340, 0.7, '6', color=colors[2], fontsize=fs, ha='left', va='top')
    axisB.plot([left+10, right, right, left+10, left+10], [bottom, bottom, top, top, bottom], '0.5')
    axisB.plot([left + 10+10], [asym_y], 'o', color=colors[2], markeredgecolor='none', markersize=3)
    axisB.text(left + 20+10, asym_y, 'Asymmetric, \n fixed repair', va='center', ha='left', fontsize=fs)
    axisB.plot([left + 10+10], [sym_y], 'o', color=colors[3], markeredgecolor='none', markersize=3)
    axisB.text(left + 20+10, sym_y, 'Symmetric, \n fixed repair', va='center', ha='left', fontsize=fs)
    axisC.text(290, 0.1, 'Generation 0', color=colors[4], fontsize=fs, ha='left', va='top')
    axisC.text(290, 0.32, '1', color=colors[4], fontsize=fs, ha='left', va='top')
    axisC.text(290, 0.6, '2+', color=colors[4], fontsize=fs, ha='left', va='top')
    #axisC.text(390, 0.50, '3+', color=colors[4], fontsize=fs, ha='left', va='top')
    axisC.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], '0.5')
    axisC.plot([left + 10], [asym_y], 'o', color=colors[4], markeredgecolor='none', markersize=3)
    axisC.text(left + 20, asym_y, 'Asymmetric, \n adaptive repair', va='center', ha='left', fontsize=fs)
    axisC.plot([left + 10], [sym_y], 'o', color=colors[5], markeredgecolor='none', markersize=3)
    axisC.text(left + 20, sym_y, 'Symmetric, \n adaptive repair', va='center', ha='left', fontsize=fs)
    
    os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/following_old_cell/')
    fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.1,wspace=0.1, hspace=0.1)
    fig.savefig('Fig_S4_'+str(a+1)+stype, bbox_inches='tight', dpi=600)
    plt.close()
  
times, generations, spec_growths, ages, inact_reps, act_reps, inact_gros, act_gros = [], [], [], [], [], [], [], []
for i in range(18):
    name = 'old_pole/'+names[i]
    time, generation, spec_growth, age, inact_rep, act_rep, inact_gro, act_gro = build_life_history(name)
    lists = [generation, spec_growth, age, inact_rep, act_rep, inact_gro, act_gro]
    time1 = time
    for j in range(7):
        list1, list2 = time1, lists[j]
        time, lists[j] = (list(x) for x in zip(*sorted(zip(list1, list2), key=itemgetter(0))))
    for k in range(len(time1)):
        time[k] = float(time[k])
        lists[0][k] = float(lists[0][k])
        lists[1][k] = float(lists[1][k])
        lists[2][k] = float(lists[2][k])
        lists[3][k] = float(lists[3][k])
        lists[4][k] = float(lists[4][k])
        lists[5][k] = float(lists[5][k])
        lists[6][k] = float(lists[6][k])
    times.append(time)
    generations.append(lists[0])
    spec_growths.append(lists[1])
    ages.append(lists[2])
    inact_reps.append(lists[3])
    act_reps.append(lists[4])
    inact_gros.append(lists[5])
    act_gros.append(lists[6])
all_bins, all_frequencies, spec_growths_100, spec_growths_100_all, spec_growth_error = [], [], [], [], []

for j in range(9):
    simulation = 'steady_state/'+names2[j]
    bins, frequencies, spec_growth = build_population_structure(simulation)
    am = numpy.mean(spec_growth) #gmean
    amstd = numpy.std(spec_growth)
    spec_growths_100.append(am)
    spec_growths_100_all.append(spec_growth)
    spec_growth_error.append(amstd)
    all_bins.append(bins)
    all_frequencies.append(frequencies)

def final_plot():
    xlim, ylim = 26, 0.6
    #AS_AR, AS_OB, AS_NR, S_AR, S_OB, S_NR
    #new = AS_AR, AS_NR, AS_OB, S_AR, S_NR, S_OB
    colorsAA = ['#fbbf07', 'red', '#F781F3', 'green', 'blue', '#81BEF7']
    colorsBB = ['#fbbf07', 'red', '#F781F3']
    colorsA = ['#fbbf07', '#F781F3', 'red', 'green', '#81BEF7', 'blue']
    colorsB = ['#fbbf07', '#F781F3', 'red']
    fs, fs2, ls = 10, 10, 1.2
    for b in range(3):
        fig = plt.figure(figsize=(8.27, 4))
        axisA, axisB, cax = fig.add_subplot(131), fig.add_subplot(132), fig.add_subplot(133)
        axisA.set_title('A', loc='left'), axisB.set_title('B', loc='left') 
        pad = 0.05
        side = "right"
        cax.set_xticklabels(['']*10)
        cax.set_yticklabels(['']*10)
        for spine in ['right', 'left', 'top', 'bottom']:
            cax.spines[spine].set_color('none')
        cax.tick_params(top="off", bottom="off", right="off", left="off")
        x = [0.05, 0.2]
        asnry = [0.9, 0.9]
        cax.plot(x, asnry, color=colorsA[2])
        cax.text(0.25, 0.89, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.845, 'no repair', fontsize=fs2)
        asory = [0.78, 0.78]
        cax.plot(x, asory, color=colorsA[1])
        cax.text(0.25, 0.77, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.725, 'fixed repair', fontsize=fs2)
        asdry = [0.66, 0.66]
        cax.plot(x, asdry, color=colorsA[0])
        cax.text(0.25, 0.65, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.61, 'adaptive repair', fontsize=fs2)
        snry = [0.545, 0.545]
        cax.plot(x, snry, color=colorsA[5])
        cax.text(0.25, 0.535, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.49, 'no repair', fontsize=fs2)
        sory = [0.425, 0.425]
        cax.plot(x, sory, color=colorsA[4])
        cax.text(0.25, 0.415, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.37, 'fixed repair', fontsize=fs2)
        sdry = [0.305, 0.305]
        cax.plot(x, sdry, color=colorsA[3])
        cax.text(0.25, 0.295, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.25, 'adaptive repair', fontsize=fs2)
        cax.set_xlim([0,1])
        cax.set_ylim([0,1])
        for a in range(6):
            number = (a*3)+b
            time = times[number]
            growth = spec_growths[number]
            generation = generations[number]
            x, y = [], []
            for c in range(len(time)):
                if time[c] != 0:
                    x.append(time[c])
                    y.append(growth[c])
                #print growth[c+1], growth[c]
                if c == len(time)-1:
                    axisA.plot(x, y, color=colorsAA[a], linewidth=ls)
                elif a < 3:
                    if generation[c+1] > generation[c]:
                        axisA.plot(x, y, color=colorsAA[a], linewidth=ls)
                        x, y = [], []
        xticks = [0, 4, 8, 12, 16, 20, 24]
        xticks_72 = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]
        axisA.set_xticks(xticks)
        axisA.set_xlabel('Time (h)', fontsize=fs)
        axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)', fontsize=fs)
        axisA.set_title('Following Old Pole Cell', fontsize=fs)
        axisA.set_xlim([0, xlim])
        axisA.set_ylim([0, ylim])
        if b == 0:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(2, 0.31, '1', color='#663300')
            axisA.text(2, 0.2, '2', color='#663300')
            axisA.text(8, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(10, 0.07, '3', color=colorsA[1]) #OB
            axisA.text(5.5, 0.01, '3', color=colorsA[2]) #NR
            axisA.text(12, 0.19, '4', color=colorsA[0])
            axisA.text(12, 0.02, '4', color=colorsA[1]) #OB
            axisA.text(15.5, 0.21, '5', color=colorsA[0])
            axisA.text(19, 0.21, '6', color=colorsA[0])
            axisA.text(22, 0.21, '7', color=colorsA[0])
            axisA.text(25, 0.21, '8', color=colorsA[0])
            axisB.text(0.04, 0.5, '0', color='#663300')
            axisB.text(0.07, 0.32, '1', color='#663300')
            axisB.text(0.025, 0.25, '2', color=colorsA[1])
            axisB.text(0.012, 0.12, '3', color=colorsA[1])
            axisB.text(0.03, 0.19, '2+', color=colorsA[0])
            axisB.text(0.037, 0.14, '2', color=colorsA[2])
        elif b == 1:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(2, 0.28, '1', color='#663300')
            axisA.text(2.5, 0.17, '2', color='#663300')
            axisA.text(10, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(12, 0.07, '3', color=colorsA[1]) #OB
            axisA.text(6.5, 0.08, '2', color=colorsA[2]) #NR
            axisA.text(13, 0.19, '4', color=colorsA[0])
            axisA.text(16.5, 0.21, '5', color=colorsA[0])
            axisA.text(20, 0.21, '6', color=colorsA[0])
            axisA.text(23, 0.21, '7', color=colorsA[0])
            #axisA.text(26, 0.21, '8', color=colorsA[0])
            axisB.text(0.04, 0.5, '0', color='#663300')
            axisB.text(0.07, 0.32, '1', color='#663300')
            axisB.text(0.025, 0.23, '2', color=colorsA[1])
            axisB.text(0.012, 0.12, '3', color=colorsA[1])
            axisB.text(0.028, 0.17, '2+', color=colorsA[0])
            axisB.text(0.035, 0.14, '2', color=colorsA[2])
        elif b == 2:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(1.5, 0.33, '1', color='#663300')
            axisA.text(4, 0.25, '2', color='#663300')
            axisA.text(3.5, 0.13, '2', color=colorsA[2])
            axisA.text(8, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(7, 0.08, '3', color=colorsA[1]) #OB
            axisA.text(5.5, 0.01, '3', color=colorsA[2]) #NR
            axisA.text(11, 0.19, '4', color=colorsA[0])
            axisA.text(18, 0.04, '4', color=colorsA[1]) #OB
            axisA.text(15, 0.21, '5', color=colorsA[0])
            axisA.text(18, 0.21, '6', color=colorsA[0])
            axisA.text(22, 0.21, '7', color=colorsA[0])
            axisA.text(25, 0.21, '8', color=colorsA[0])
            axisB.text(0.035, 0.5, '0', color='#663300')
            axisB.text(0.06, 0.32, '1', color='#663300')
            axisB.text(0.027, 0.25, '2', color=colorsA[1])
            axisB.text(0.013, 0.12, '3', color=colorsA[1])
            axisB.text(0.025, 0.17, '2+', color=colorsA[0])
            axisB.text(0.04, 0.14, '2', color=colorsA[2])
            axisB.text(0.04, 0.02, '3', color=colorsA[2])
        plt.setp( axisB.get_yticklabels(), visible=False)
        axisB.set_xlim([0, 0.105])
        axisB.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10])
        axisB.set_ylim([0, 0.6])
        axisB.set_xlabel('Frequency in population', fontsize=fs)
        axisB.set_title('Growth Rate Distribution', fontsize=fs)
        for d in range(3):
            number = (d*3)+b
            bins, frequency = all_bins[number], all_frequencies[number]
            axisB.plot(frequency, bins, color=colorsBB[d])
        pad = 0.05
        side = "left"
        divider = make_axes_locatable(cax)
        caxB = divider.append_axes(side, size="25%", pad=0, axisbg='none', frameon=True)
        caxB.set_xticklabels(['']*10)
        caxB.set_yticklabels(['']*10)
        #caxB.spines['left'].set_color('none')
        caxB.tick_params(top="off", right="off", bottom="off", left="off")
        x = [1.5, 1, 0.5]
        for_plotting = []
        for e in range(3):
            number = (e*3)+b
            mean = spec_growths_100[number]
            std = spec_growth_error[number]
            all_spec = spec_growths_100_all[number]
            for_plotting.append(all_spec)
            #caxB.errorbar(x[e], mean, yerr=std, marker='o', color=colorsB[e])
        whiskerprops = dict(linestyle='-.', linewidth=1, color='k')
        medianprops = dict(linestyle='-', linewidth=2, color='k')
        bp = caxB.boxplot(for_plotting, positions=[1, 2, 3], medianprops=medianprops, whiskerprops=whiskerprops)
        for box in range(len(bp['boxes'])):
            bp['boxes'][box].set(color=colorsBB[box], linewidth=2)
        caxB.set_xlim([0, 4])
        caxB.set_ylim([0, 0.6])
        #caxB.set_title(r'Means $\pm$ SD', fontsize=fs)
        caxB.set_title('Boxplots', fontsize=fs)
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/following_old_cell/')
        fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.05, hspace=0.2)
        fig.savefig('Fig_2_'+str(b+1)+stype, bbox_inches='tight', dpi=600)
        plt.close()
        
        fig = plt.figure(figsize=(8.27, 4))
        axisC = fig.add_subplot(111)
        side = "right"
        divider = make_axes_locatable(axisC)
        cax = divider.append_axes(side, size="25%", pad=pad, axisbg='none')
        cax.set_xticklabels(['']*10)
        cax.set_yticklabels(['']*10)
        for spine in ['right', 'left', 'top', 'bottom']:
            cax.spines[spine].set_color('none')
        cax.tick_params(top="off", bottom="off", right="off", left="off")
        x = [0.05, 0.2]
        asnry = [0.9, 0.9]
        cax.plot(x, asnry, color=colorsA[0])
        cax.text(0.25, 0.89, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.845, 'adaptive repair', fontsize=fs2)
        asory = [0.78, 0.78]
        cax.plot(x, asory, color=colorsA[3])
        cax.text(0.25, 0.77, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.725, 'adaptive repair', fontsize=fs2)
        cax.set_xlim([0,1])
        cax.set_ylim([0,1])
        axisC.set_ylabel(r'$ \^\beta $', style='italic', fontsize=fs)
        axisC.set_xticks(xticks_72)
        axisC.set_xlabel('Time (h)', fontsize=fs)
        axisC.set_title('Following old pole cell over time', fontsize=fs)
        axisC.set_xlim([0, 72])
        colorsC = [colorsA[0], colorsA[3], colorsA[4]]
        muS = 0.6
        Yr = 0.8
        numbers = [[0, 1, 2], [9, 10, 11], [12, 13, 14]]
        for i in range(2):
            number = numbers[i][b]
            time, age_using, generation = times[number], ages[number], generations[number]
            x, y = [], []
            for j in range(len(time)):
                x.append(time[j])
                age = age_using[j]
                if age == 0:
                    y.append(0)
                elif age >= 1:
                    y.append(1)
                else:
                    age_appending = float(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
                    y.append(age_appending)
                    if j == len(time)-1:
                        axisC.plot(x, y, color=colorsC[i], linewidth=ls)
                    elif i == 1:
                        continue
                    elif generation[j+1] > generation[j]:
                        axisC.plot(x, y, color=colorsC[i], linewidth=ls)
                        x, y = [], []
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/following_old_cell/')
        fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.0, hspace=0.2)
        fig.savefig('Fig_1_'+str(b+1)+stype, bbox_inches='tight', dpi=600)
        plt.close()
final_plot()

def combined_plot():
    xlim, ylim = 26, 0.6
    #AS_AR, AS_OB, AS_NR, S_AR, S_OB, S_NR
    #new = AS_AR, AS_NR, AS_OB, S_AR, S_NR, S_OB
    colorsAA = ['#fbbf07', 'red', '#F781F3', 'green', 'blue', '#81BEF7']
    colorsBB = ['#fbbf07', 'red', '#F781F3']
    colorsA = ['#fbbf07', '#F781F3', 'red', 'green', '#81BEF7', 'blue']
    colorsB = ['#fbbf07', '#F781F3', 'red']
    fs, fs2, ls = 10, 10, 1.2
    for b in range(3):
        fig = plt.figure(figsize=(8.27, 8.27))
        axisA, axisB = fig.add_subplot(223), fig.add_subplot(224)
        axisA.set_title('B', loc='left'), axisB.set_title('C', loc='left') 
        
        for a in range(6):
            number = (a*3)+b
            time = times[number]
            growth = spec_growths[number]
            generation = generations[number]
            x, y = [], []
            for c in range(len(time)):
                if time[c] != 0:
                    x.append(time[c])
                    y.append(growth[c])
                #print growth[c+1], growth[c]
                if c == len(time)-1:
                    axisA.plot(x, y, color=colorsAA[a], linewidth=ls)
                elif a < 3:
                    if generation[c+1] > generation[c]:
                        axisA.plot(x, y, color=colorsAA[a], linewidth=ls)
                        x, y = [], []
        xticks = [0, 4, 8, 12, 16, 20, 24]
        xticks_72 = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]
        axisA.set_xticks(xticks)
        axisA.set_xlabel('Time (h)', fontsize=fs)
        axisA.set_ylabel(r'Cellular specific growth rate (h$^{-1}$)', fontsize=fs)
        axisA.set_title('Following Old Pole Cell', fontsize=fs)
        axisA.set_xlim([0, xlim])
        axisA.set_ylim([0, ylim])
        if b == 0:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(2, 0.31, '1', color='#663300')
            axisA.text(2, 0.2, '2', color='#663300')
            axisA.text(8, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(10, 0.07, '3', color=colorsA[1]) #OB
            axisA.text(5.5, 0.01, '3', color=colorsA[2]) #NR
            axisA.text(12, 0.19, '4', color=colorsA[0])
            axisA.text(12, 0.02, '4', color=colorsA[1]) #OB
            axisA.text(15.5, 0.21, '5', color=colorsA[0])
            axisA.text(19, 0.21, '6', color=colorsA[0])
            axisA.text(22, 0.21, '7', color=colorsA[0])
            axisA.text(25, 0.21, '8', color=colorsA[0])
            axisB.text(0.04, 0.5, '0', color='#663300')
            axisB.text(0.07, 0.32, '1', color='#663300')
            axisB.text(0.025, 0.25, '2', color=colorsA[1])
            axisB.text(0.012, 0.12, '3', color=colorsA[1])
            axisB.text(0.03, 0.19, '2+', color=colorsA[0])
            axisB.text(0.037, 0.14, '2', color=colorsA[2])
        elif b == 1:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(2, 0.28, '1', color='#663300')
            axisA.text(2.5, 0.17, '2', color='#663300')
            axisA.text(10, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(12, 0.07, '3', color=colorsA[1]) #OB
            axisA.text(6.5, 0.08, '2', color=colorsA[2]) #NR
            axisA.text(13, 0.19, '4', color=colorsA[0])
            axisA.text(16.5, 0.21, '5', color=colorsA[0])
            axisA.text(20, 0.21, '6', color=colorsA[0])
            axisA.text(23, 0.21, '7', color=colorsA[0])
            #axisA.text(26, 0.21, '8', color=colorsA[0])
            axisB.text(0.04, 0.5, '0', color='#663300')
            axisB.text(0.07, 0.32, '1', color='#663300')
            axisB.text(0.025, 0.23, '2', color=colorsA[1])
            axisB.text(0.012, 0.12, '3', color=colorsA[1])
            axisB.text(0.028, 0.17, '2+', color=colorsA[0])
            axisB.text(0.035, 0.14, '2', color=colorsA[2])
        elif b == 2:
            axisA.text(1.8, 0.5, '0', color='#663300')
            axisA.text(1.5, 0.33, '1', color='#663300')
            axisA.text(4, 0.25, '2', color='#663300')
            axisA.text(3.5, 0.13, '2', color=colorsA[2])
            axisA.text(8, 0.18, '3', color=colorsA[0]) #AR
            axisA.text(7, 0.08, '3', color=colorsA[1]) #OB
            axisA.text(5.5, 0.01, '3', color=colorsA[2]) #NR
            axisA.text(11, 0.19, '4', color=colorsA[0])
            axisA.text(18, 0.04, '4', color=colorsA[1]) #OB
            axisA.text(15, 0.21, '5', color=colorsA[0])
            axisA.text(18, 0.21, '6', color=colorsA[0])
            axisA.text(22, 0.21, '7', color=colorsA[0])
            axisA.text(25, 0.21, '8', color=colorsA[0])
            axisB.text(0.035, 0.5, '0', color='#663300')
            axisB.text(0.06, 0.32, '1', color='#663300')
            axisB.text(0.027, 0.25, '2', color=colorsA[1])
            axisB.text(0.013, 0.12, '3', color=colorsA[1])
            axisB.text(0.025, 0.17, '2+', color=colorsA[0])
            axisB.text(0.04, 0.14, '2', color=colorsA[2])
            axisB.text(0.04, 0.02, '3', color=colorsA[2])
        plt.setp( axisB.get_yticklabels(), visible=False)
        axisB.set_xlim([0, 0.105])
        axisB.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10])
        axisB.set_ylim([0, 0.6])
        axisB.set_xlabel('Frequency in population', fontsize=fs)
        axisB.set_title('Growth Rate Distribution', fontsize=fs)
        for d in range(3):
            number = (d*3)+b
            bins, frequency = all_bins[number], all_frequencies[number]
            axisB.plot(frequency, bins, color=colorsBB[d])
        pad = 0.05
        side = "right"
        divider = make_axes_locatable(axisB)
        caxB = divider.append_axes(side, size="25%", pad=0, axisbg='none', frameon=True)
        caxB.set_xticklabels(['']*10)
        caxB.set_yticklabels(['']*10)
        #caxB.spines['left'].set_color('none')
        caxB.tick_params(top="off", right="off", bottom="off", left="off")
        x = [1.5, 1, 0.5]
        for_plotting = []
        for e in range(3):
            number = (e*3)+b
            mean = spec_growths_100[number]
            std = spec_growth_error[number]
            all_spec = spec_growths_100_all[number]
            for_plotting.append(all_spec)
            #caxB.errorbar(x[e], mean, yerr=std, marker='o', color=colorsB[e])
        whiskerprops = dict(linestyle='-.', linewidth=1, color='k')
        medianprops = dict(linestyle='-', linewidth=2, color='k')
        bp = caxB.boxplot(for_plotting, positions=[1, 2, 3], medianprops=medianprops, whiskerprops=whiskerprops)
        for box in range(len(bp['boxes'])):
            bp['boxes'][box].set(color=colorsBB[box], linewidth=2)
        caxB.set_xlim([0, 4])
        caxB.set_ylim([0, 0.6])
        #caxB.set_title(r'Means $\pm$ SD', fontsize=fs)
        caxB.set_title('Boxplots', fontsize=fs)
        #os.chdir('/Users/u1560915/git/iDynoMiCS/paper_july_2018/following_old_cell/')
        #fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.05, hspace=0.2)
        #fig.savefig('Fig_R2_'+str(b+1)+stype, bbox_inches='tight')
        #plt.close()
        
        axisC = fig.add_subplot(211)
        axisC.set_title('A', loc='left')
        side = "right"
        divider = make_axes_locatable(axisC)
        cax = divider.append_axes(side, size="25%", pad=pad, axisbg='none')
        cax.set_xticklabels(['']*10)
        cax.set_yticklabels(['']*10)
        for spine in ['right', 'left', 'top', 'bottom']:
            cax.spines[spine].set_color('none')
        cax.tick_params(top="off", bottom="off", right="off", left="off")
        x = [0.05, 0.2]
        asnry = [0.9, 0.9]
        
        cax.plot(x, asnry, color=colorsA[2])
        cax.text(0.25, 0.89, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.845, 'no repair', fontsize=fs2)
        asory = [0.78, 0.78]
        cax.plot(x, asory, color=colorsA[1])
        cax.text(0.25, 0.77, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.725, 'fixed repair', fontsize=fs2)
        asdry = [0.66, 0.66]
        cax.plot(x, asdry, color=colorsA[0])
        cax.text(0.25, 0.65, 'Asymmetric,', fontsize=fs2)
        cax.text(0.26, 0.61, 'adaptive repair', fontsize=fs2)
        snry = [0.545, 0.545]
        cax.plot(x, snry, color=colorsA[5])
        cax.text(0.25, 0.535, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.49, 'no repair', fontsize=fs2)
        sory = [0.425, 0.425]
        cax.plot(x, sory, color=colorsA[4])
        cax.text(0.25, 0.415, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.37, 'fixed repair', fontsize=fs2)
        sdry = [0.305, 0.305]
        cax.plot(x, sdry, color=colorsA[3])
        cax.text(0.25, 0.295, 'Symmetric,', fontsize=fs2)
        cax.text(0.26, 0.25, 'adaptive repair', fontsize=fs2)
        cax.set_xlim([0,1])
        cax.set_ylim([0,1])
        
        axisC.set_ylabel(r'$ \^\beta $', style='italic', fontsize=fs)
        axisC.set_xticks(xticks_72)
        axisC.set_xlabel('Time (h)', fontsize=fs)
        axisC.set_title('Following old pole cell over time', fontsize=fs)
        axisC.set_xlim([0, 72])
        colorsC = [colorsA[0], colorsA[3], colorsA[4]]
        muS = 0.6
        Yr = 0.8
        numbers = [[0, 1, 2], [9, 10, 11], [12, 13, 14]]
        for i in range(2):
            number = numbers[i][b]
            time, age_using, generation = times[number], ages[number], generations[number]
            x, y = [], []
            for j in range(len(time)):
                x.append(time[j])
                age = age_using[j]
                if age == 0:
                    y.append(0)
                elif age >= 1:
                    y.append(1)
                else:
                    age_appending = float(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
                    y.append(age_appending)
                    if j == len(time)-1:
                        axisC.plot(x, y, color=colorsC[i], linewidth=ls)
                    elif i == 1:
                        continue
                    elif generation[j+1] > generation[j]:
                        axisC.plot(x, y, color=colorsC[i], linewidth=ls)
                        x, y = [], []
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/following_old_cell/')
        fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.0, hspace=0.2)
        fig.savefig('Fig_1_2_'+str(b+1)+stype, bbox_inches='tight', dpi=600)
        plt.close()
combined_plot()

def prep_ptot():
    fs = 10
    colors = ['#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#06e4ea', '#0066FF', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea', '#FF00FF', '#FF0000', '#FF9900', '#f4de06', '#00FF00', '#0066FF','#06e4ea']
    muS = 0.6
    Yr = 0.8
    for b in range(3):
        fig = plt.figure(figsize=(8.27, 8.27))
        axisA, axisB = fig.add_subplot(221), fig.add_subplot(222)
        axisA.set_title('A', loc='left'), axisB.set_title('B', loc='left')
        for a in range(2):
            if a == 1:
                a = 3
            number = (a*3)+b
            gen = generations[number]
            inactGro = inact_gros[number]
            inactRep = inact_reps[number]
            actGro = act_gros[number]
            actRep = act_reps[number]
            age = ages[number]
            linex, liney  = [0, 0.7], [0, 0.7]
            x, y = [], []
            oldx, oldy = [], []
            current_gen = 0
            if a == 0:
                axisone = axisA
                axisone.set_xlim([0, 0.7])
                axisone.set_ylim([0, 0.7])
                axisone.set_title('Asymmetric')
                axisone.set_xlabel(r'$ \^\beta $', style='italic', fontsize=fs)
                axisone.set_ylabel(r'$P_{rep}/P_{tot}$', fontsize=fs)
            else:
                axisone = axisB
                axisone.set_xlim([0, 0.1])
                axisone.set_ylim([0, 0.1])
                axisone.set_title('Symmetric')
                axisone.set_xlabel(r'$ \^\beta $', style='italic', fontsize=fs)
            for c in range(0, len(age)):
                pActGro = actGro[c]
                pActRep = actRep[c]
                pInactGro = inactGro[c]
                pInactRep = inactRep[c]
                pRep = (pActRep + pInactRep)
                pTot = (pActGro + pActRep + pInactGro + pInactRep)
                axisone.plot(linex, liney, 'k')
                if gen[c] > current_gen:
                    axisone.plot(x, y, colors[current_gen])
                    oldx.append(x[-1])
                    oldy.append(y[-1])
                    x, y = [], []
                if current_gen < 7:
                    x.append([float(( age[c] / (1 - age[c]))  * ( math.sqrt((Yr / muS) * (1 / (1 - age[c]))) - 1))])
                    y.append([pRep/pTot])
                    if gen[c] > current_gen:
                        w = x[0]
                        wx = ([w, oldx[0]])
                        z = y[0]
                        zy = ([z, oldy[0]])
                        axisone.plot(wx, zy, linestyle = '--', color = '#565353')
                        axisone.plot(oldx[0], oldy[0], marker='o', color=colors[current_gen])
                        current_gen += 1
                        oldx, oldy= [], []
                else:
                    break
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/following_old_cell/')
        fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.2, hspace=0.2)
        fig.savefig('Fig_S1_'+str(b+1)+stype, bbox_inches='tight', dpi=600)
        plt.close()
    return
prep_ptot()
