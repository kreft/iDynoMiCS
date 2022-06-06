#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting_age as toolbox_plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable
import toolbox_results
import matplotlib
import numpy
import math
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import calc_roughness
import matplotlib as mpl
from datetime import datetime
startTime = datetime.now()

treatments = ['SAR_ASNR', 'ASNR', 'SAR']
fs, fs2, wspace, hspace, lb = 10, 10, 0.1, 0.4, 0.06
xticks = [0, 50, 100, 150, 200, 250]

def get_highest_growth_rate(agent_output1, agent_output2, agent_output3):
    agent_outputs = [agent_output1, agent_output2, agent_output3]
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    all_growth = []
    for a in range(3):
        growthvec = []
        for species in agent_outputs[a].species_outputs:
            for cell in species.members:
                growthrate = cell.get_specific_growth_rate(biomass_names)
                growthrate = float(growthrate)
                growthvec.append(growthrate)
            all_growth.append(max(growthvec))
    return max(all_growth)

def color_cells_growth_rate(agent_output):
    species_color_dict = {'OldieA' : 'jet_r'}
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    print ('Coloring cells')
    time = (agent_output.time)
    for species in agent_output.species_outputs:
        if species.members == []:
            continue
        for cell in species.members:
            name = species.name
            for (species_name, colormap_name) in species_color_dict.iteritems():
                if name == species_name:
                    norm = matplotlib.colors.Normalize(vmin=0,vmax=0.6)
                    colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                    growth = float(cell.get_specific_growth_rate(biomass_names))
                cell.color=m.to_rgba(growth)
    return time
    
def plot(iter_info, min_max_concns, b, c, axis, roughness, label):
    axis = axis[b,c]
    time = color_cells_growth_rate(iter_info.agent_output)
    toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
    axis.fill_between([0, nJ*res], [0]*2, y2=[-res]*2, color='k', zorder=-1)
    if min_max_concns != {}:
        solute_output = toolbox_results.SoluteOutput(iter_info.env_output, name="glucose")
        toolbox_idynomics.solute_contour(axis, solute_output, concn_range=min_max_concns["glucose"],interpolation='bicubic')
    if b == 2:
        axis.set_xlabel(r'x ($\mu$m)', fontsize=fs)
    if c == 0:
        axis.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    axis.set_xticklabels(xticks, visible=True, fontsize=8)
    axis.set_yticklabels(xticks, visible=True, fontsize=8)
    axis.text(250, ((nI+0.1)*res-8)-20, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=8)
    axis.text(250, ((nI+0.1)*res-8)-20, r'$\sigma_{f}$ = %.1f'%float(roughness), va='bottom', ha='right', color='#bd0303', fontsize=8)
    axis.text(30, ((nI+0.1)*res-8)-20, label, va='bottom', ha='right', color='k', fontsize=10)
    axis.set_xticks(xticks)
    axis.set_xlim(0, nJ * res)
    axis.set_ylim(0, nI * res)
    return
    
os.chdir('/Volumes/Robyn_W_2/july_2018/paper/roughness/')

low = ['NANRL152_25day(20160315_1337)', 'NANRL152_25day(20160317_1534)', 'NANRL152_25day(20160319_1939)', 'NANRL152_25day(20160321_1746)']
med = ['NANRM152_25day(20160316_1440)', 'NANRM152_25day(20160318_1835)', 'NANRM152_25day(20160320_1515)', 'NANRM152_25day(20160322_1524)']
high = ['NANRH152_25day(20160315_0935)', 'NANRH152_25day(20160317_1147)', 'NANRH152_25day(20160319_1645)', 'NANRH152_25day(20160321_1108)']
biofilms = [high, med, low]

def colorbar_figure():
    fig = plt.figure(figsize=(8.27, 4))
    ax1 = fig.add_subplot(8,2,1)
    ax2 = fig.add_subplot(8,2,3)
    cmap = mpl.cm.jet_r
    norm = mpl.colors.Normalize(vmin=0, vmax=0.6)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
    cb1.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    cb1.ax.tick_params(labelsize=8)
    cb1.ax.text(0.5, 0.5, 'Specific growth rate', ha='center', va='center', fontsize=8)
    cmap = mpl.cm.gray
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal')
    cb2.set_ticks([0, 1])
    cb2.set_ticklabels(['No solute', 'Concentrated \n solute'])
    cb2.ax.tick_params(labelsize=8)
    plt.tight_layout()
    fig.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.5, hspace=1.5)
    os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/roughness/')
    plt.savefig('colorbars.pdf', dpi = 200)
    return


for a in range(1): #4
    os.chdir('/Volumes/Robyn_W_2/july_2018/paper/roughness/')
    for b in range(1): #3
        file_name = biofilms[a][b]
        sim = toolbox_idynomics.SimulationDirectory(file_name)
        last_iter = sim.get_last_iterate_number()
        last_iter -= (last_iter % 20)
        height, roughness = [], []
        for c in range(int(last_iter/20)):
            num = c*20
            r, h = calc_roughness.calc_roughness(file_name, num)
            roughness.append(r)
            height.append(h*4)
        

    
    

for a in range(4):
    #a += 1
    os.chdir('/Volumes/Robyn_W_2/july_2018/paper/roughness/')
    f, axis = plt.subplots(4,3, sharex=True, sharey=True) #r, c
    labels = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I']]
    for b in range(3):
        file_name = biofilms[b][a]
        sim = toolbox_idynomics.SimulationDirectory(file_name)
        num_digits = len(str(sim.get_last_iterate_number()))
        nI, nJ, nK, res = sim.find_domain_dimensions()        
        last_iter = sim.get_last_iterate_number()
        first, second, third = last_iter/3, (last_iter/3)*2, last_iter
        first -= (first % 20)
        second -= (second % 20)
        iter_info_1, iter_info_2, iter_info_3 = sim.get_single_iterate(first), sim.get_single_iterate(second), sim.get_single_iterate(third)
        min_max_concns_1, min_max_concns_2, min_max_concns_3 = iter_info_1.get_min_max_concns(), iter_info_2.get_min_max_concns(), iter_info_3.get_min_max_concns()
        print min_max_concns_1
        print min_max_concns_2
        print min_max_concns_3
        #max_growth = get_highest_growth_rate(iter_info_1.agent_output, iter_info_2.agent_output, iter_info_3.agent_output)
        roughness_1, roughness_2, roughness_3 = calc_roughness.calc_roughness(file_name, first)[0], calc_roughness.calc_roughness(file_name, second)[0], calc_roughness.calc_roughness(file_name, third)[0]
        plot(iter_info_1, min_max_concns_1, b, 0, axis, roughness_1, labels[b][0])
        plot(iter_info_2, min_max_concns_2, b, 1, axis, roughness_2, labels[b][1])
        plot(iter_info_3, min_max_concns_3, b, 2, axis, roughness_3, labels[b][2])
        sim.clean_up()
    axis[3,0].axis('off')
    axis[3,1].axis('off')
    axis[3,2].axis('off')
    
    divider = make_axes_locatable(axis[3,1])
    
    caxB = divider.append_axes("top", size="20%", pad=0.1, axisbg='none', frameon=False)
    caxB.set_xticklabels(['']*10)
    caxB.set_yticklabels(['']*10)
    caxB.spines['left'].set_color('none')
    caxB.tick_params(top="off", right="off", bottom="off", left="off")
    cmap = mpl.cm.gray
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb2 = mpl.colorbar.ColorbarBase(caxB, cmap=cmap, norm=norm, orientation='horizontal')
    cb2.set_ticks([0, 1])
    cb2.set_ticklabels(['No solute', 'Concentrated \n solute'])
    cb2.ax.tick_params(labelsize=8)
    
    caxC = divider.append_axes("top", size="20%", pad=0.25, axisbg='none', frameon=False)
    caxC.set_xticklabels(['']*10)
    caxC.set_yticklabels(['']*10)
    caxC.spines['left'].set_color('none')
    caxC.tick_params(top="off", right="off", bottom="off", left="off")
    cmap = mpl.cm.jet_r
    norm = mpl.colors.Normalize(vmin=0, vmax=0.6)
    cb1 = mpl.colorbar.ColorbarBase(caxC, cmap=cmap, norm=norm, orientation='horizontal')
    cb1.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    cb1.ax.tick_params(labelsize=8)
    cb1.ax.text(0.5, 0.5, 'Specific growth rate', ha='center', va='center', fontsize=8)    
    
    axis[2,1].set_xlabel(r'x ($\mu$m)'+'\n ', fontsize=fs)
    axis[2,1].set_xticklabels(xticks, visible=True, fontsize=8)
    axis[2,1].set_xticks(xticks)
    f.set_size_inches(8.27,11)
    f.subplots_adjust(hspace=0.3)
    plt.gca()
    #.set_aspect('equal')
    os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/roughness/')
    plt.savefig('Testing_colorbar.png', dpi=200)
    plt.savefig('roughness_'+str(a+1)+'.png', dpi = 200)
    plt.close()
    print(datetime.now() - startTime)
print(datetime.now() - startTime)

