#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 12:17:01 2020

@author: robynwright
"""

from __future__ import division
from __future__ import with_statement
import os
import toolbox_basic as basic
import toolbox_idynomics
import toolbox_results as results
from mpl_toolkits.axes_grid1 import make_axes_locatable
import toolbox_results
import matplotlib
import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
from operator import itemgetter
from matplotlib.patches import Patch

fs, fs2, wspace, hspace, lb = 10, 12, 0.1, 0.4, 0.06
xticks = [0, 50, 100, 150, 200, 250]

def plot_biofilm(ax, fi, arate=True):
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair','inactiveBiomassGrowth', 'inactiveBiomassRepair']
    species_color_dict = {'OldieA' : 'cool', 'OldieB' : 'autumn'}
    os.chdir(fi+'/lastIter')
    output = results.AgentOutput(path='agent_State(last).xml')
    time = (output.time)
    ma1, ma2 = 0, 0
    for species in output.species_outputs:
        if species.members == []:
            continue
        if not arate:
            all_growth = []
            for cell in species.members:
                all_growth.append(float(cell.get_specific_growth_rate(biomass_names)))
            ma = max(all_growth)
            if ma1 == 0: ma1 = ma
            else: ma2 = ma
        for cell in species.members:
            name = species.name
            for (species_name, colormap_name) in species_color_dict.items():
                if name == species_name:
                    colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                    if not arate:
                        l0, l06 = numpy.log(0.001), numpy.log(0.6)
                        norm = matplotlib.colors.Normalize(vmin=l0,vmax=l06)
                    else:
                        norm = matplotlib.colors.Normalize(vmin=0,vmax=1)
                    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                    if arate:
                        age = float(cell.vars['age'])
                    else:
                        gr = float(cell.get_specific_growth_rate(biomass_names))
                        if gr == 0: gr = numpy.log(0.001)
                        else: gr = numpy.log(gr)
                if arate:
                    cell.color=m.to_rgba(age)
                else:
                    cell.color=m.to_rgba(gr)
    toolbox_idynomics.plot_cells_2d(ax, output)
    output = results.EnvOutput(path='env_State(last).xml')
    for solute in output.solute_outputs:
        if solute.name == 'glucose':
            solute_output = solute
            toolbox_idynomics.solute_contour(ax, solute_output, interpolation='bicubic')
    ax.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    print(ma1, ma2)
    return

def get_colorbar(ax, color, name, side, arate=True):
    pad=0.01
    plt.sca(ax)
    plt.xticks([]), plt.yticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size="80%", pad=pad, frameon=True)
    if name != 'Solute':
        cmap = color
        if arate:
            cax.text(-0.3, 0.5, 'Age '+name, va='center', ha='center', fontsize=fs, rotation=90, transform=cax.transAxes)
            cax.text(1.5, 0.1, 'Young cells', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
            cax.text(1.5, 0.9, 'Old cells', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
        else:
            cax.text(-0.3, 0.5, name+r' specific growth rate (h$^{-1}$)', va='center', ha='center', fontsize=fs, rotation=90, transform=cax.transAxes)
            cax.text(1.5, 0.1, '<0.001', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
            cax.text(1.5, 0.5, '0.025', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
            cax.text(1.5, 0.9, '0.6+', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
    else:
        cmap = mpl.cm.gray
        cax.text(-0.3, 0.5, 'Solute concentration', va='center', ha='center', fontsize=fs, rotation=90, transform=cax.transAxes)
        cax.text(1.5, 0.1, 'No solute', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
        cax.text(1.5, 0.9, 'Concentrated\nsolute', va='center', ha='center', fontsize=fs-1, rotation=90, transform=cax.transAxes)
    if not arate:
        l0, l06 = numpy.log(0.001), numpy.log(0.6)
        norm = matplotlib.colors.Normalize(vmin=l0,vmax=l06)
    else:
        norm = matplotlib.colors.Normalize(vmin=0,vmax=1)
    cmap = matplotlib.cm.get_cmap(color, 256)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_ticks([])
    return

def get_specgrowth_ratios(folder, f, calc='mean'):
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair','inactiveBiomassGrowth', 'inactiveBiomassRepair']
    sim_dir_path = basic.check_path(folder+f+'/')
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    t_a, specgrowth_a = [], []
    t_b, specgrowth_b = [], []
    for filename in file_list:
        growth_a, growth_b = [], []
        output = results.AgentOutput(path=filename)
        t_a.append(float(output.time))
        t_b.append(float(output.time))
        for species in output.species_outputs:
            if species.members == []:
                continue
            for cell in species.members:
                name = species.name
                gr = float(cell.get_specific_growth_rate(biomass_names))
                if name == 'OldieA': growth_a.append(gr)
                elif name == 'OldieB': growth_b.append(gr)
        if calc == 'mean':
            specgrowth_a.append(numpy.mean(growth_a))
            specgrowth_b.append(numpy.mean(growth_b))
        elif calc == 'median':
            specgrowth_a.append(numpy.median(growth_a))
            specgrowth_b.append(numpy.median(growth_b))
        elif calc == 'max':
            specgrowth_a.append(max(growth_a))
            specgrowth_b.append(max(growth_b))
        elif calc == 'min':
            specgrowth_a.append(min(growth_a))
            specgrowth_b.append(min(growth_b))
    #basic.rm_dir(file_dir)
    t_a1, t_b1 = t_a, t_b
    t_a, specgrowth_a = (list(x) for x in zip(*sorted(zip(t_a1, specgrowth_a), key=itemgetter(0))))
    t_b, specgrowth_b = (list(x) for x in zip(*sorted(zip(t_b1, specgrowth_b), key=itemgetter(0))))
    specgrowth_ratio = []
    for j in range(len(t_a)):
        if float(specgrowth_b[j]) < 0.001: specgrowth_b[j] = 0.001
        if float(specgrowth_a[j]) < 0.001: specgrowth_a[j] = 0.001
        specgrowth_ratio.append(numpy.log(specgrowth_b[j]/specgrowth_a[j]))
    return t_a, specgrowth_ratio,
    

def get_biomass_ratios(folder, f):
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair','inactiveBiomassGrowth', 'inactiveBiomassRepair']
    sim_dir_path = basic.check_path(folder+f+'/')
    file_dir = os.path.join(sim_dir_path, 'agent_Sum')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    t_a, population_a, biomass_a, growthrate_a = [], [], [], []
    t_b, population_b, biomass_b, growthrate_b = [], [], [], []
    last_pop_a, last_biomass_a, last_growthrate_a = [], [], []
    last_pop_b, last_biomass_b, last_growthrate_b = [], [], []
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output, 'OldieA')
        requirements={}
        cells = species.find_cells(requirements)
        #single_result = results.SingleResult()
        t_a.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population_a.append(float(cell.vars['population']))
        biomass_a.append(float(cell.vars['mass']))
        growthrate_a.append(float(cell.vars['growthRate']))
        species = results.SpeciesOutput(output, 'OldieB')
        requirements={}
        cells = species.find_cells(requirements)
        #single_result = results.SingleResult()
        t_b.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population_b.append(float(cell.vars['population']))
        biomass_b.append(float(cell.vars['mass']))
        growthrate_b.append(float(cell.vars['growthRate']))
    basic.rm_dir(file_dir)
    lists_a, lists_b = [population_a, biomass_a, growthrate_a], [population_b, biomass_b, growthrate_b]
    t_a1, t_b1 = t_a, t_b
    for i in range(3):
        list1, list2, list3, list4 = t_a1, lists_a[i], t_b1, lists_b[i]
        t_a, lists_a[i] = (list(x) for x in zip(*sorted(zip(list1, list2), key=itemgetter(0))))
        t_b, lists_b[i] = (list(x) for x in zip(*sorted(zip(list3, list4), key=itemgetter(0))))
    biomass_ratio, population_ratio, growthrate_ratio = [], [], []
    for j in range(len(lists_a[1])):
        biomass_ratio.append(numpy.log(lists_b[1][j]/lists_a[1][j]))
        population_ratio.append(numpy.log(lists_a[0][j]/lists_b[0][j]))
        if lists_a[2][j] == 0 and lists_b[2][j] == 0:
            growthrate_ratio.append(0)
        elif lists_a[2][j] == 0:
            lists_b[2][j] = abs(1/(lists_b[2][j]))
            growthrate_ratio.append(numpy.log(lists_b[2][j]))
        elif lists_b[2][j] == 0:
            lists_a[2][j] = abs((lists_a[2][j])/1)
            growthrate_ratio.append(numpy.log(lists_a[2][j]))
        else:
            growthrate_ratio.append(abs(numpy.log(lists_a[2][j]/lists_b[2][j])))
        if j == len(lists_a[1])-1:
            last_pop_a.append(lists_a[0][j])
            last_pop_b.append(lists_b[0][j])
            last_biomass_a.append(lists_a[1][j])
            last_biomass_b.append(lists_b[1][j])
            last_growthrate_a.append(lists_a[2][j])
            last_growthrate_b.append(lists_b[2][j])
    return t_a, [biomass_ratio, population_ratio, growthrate_ratio]

def get_prop_dam_rep(folder, f):
    sim_dir_path = basic.check_path(folder+f+'/')
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    t_a, prot_oa, dam_oa, rep_oa, age_oa = [], [], [], [], []
    t_b, prot_ob, dam_ob, rep_ob, age_ob = [], [], [], [], []
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output, 'OldieA')
        requirements={}
        cells = species.find_cells(requirements)
        t_a.append(float(output.time))
        if len(cells) == 0:
            continue
        this_gro, this_rep, this_dam_gro, this_dam_rep, this_age = [], [], [], [], []
        for cell in cells:
            this_gro.append(float(cell.vars['activeBiomassGrowth']))
            this_rep.append(float(cell.vars['activeBiomassRepair']))
            this_dam_gro.append(float(cell.vars['inactiveBiomassGrowth']))
            this_dam_rep.append(float(cell.vars['inactiveBiomassRepair']))
            this_age.append(float(cell.vars['age']))
        this_gro = numpy.mean(this_gro)
        this_rep = numpy.mean(this_rep)
        this_dam_gro = numpy.mean(this_dam_gro)
        this_dam_rep = numpy.mean(this_dam_rep)
        this_age = numpy.mean(this_age)
        this_prot = this_gro+this_rep+this_dam_gro+this_dam_rep
        this_prop_dam = (this_dam_gro+this_dam_rep)/this_prot
        this_prop_rep = this_rep/this_prot
        prot_oa.append(this_prot)
        dam_oa.append(this_prop_dam)
        rep_oa.append(this_prop_rep)
        age_oa.append(this_age)
        
        species = results.SpeciesOutput(output, 'OldieB')
        requirements={}
        cells = species.find_cells(requirements)
        t_b.append(float(output.time))
        if len(cells) == 0:
            continue
        this_gro, this_rep, this_dam_gro, this_dam_rep, this_age = [], [], [], [], []
        for cell in cells:
            this_gro.append(float(cell.vars['activeBiomassGrowth']))
            this_rep.append(float(cell.vars['activeBiomassRepair']))
            this_dam_gro.append(float(cell.vars['inactiveBiomassGrowth']))
            this_dam_rep.append(float(cell.vars['inactiveBiomassRepair']))
            this_age.append(float(cell.vars['age']))
        this_gro = numpy.mean(this_gro)
        this_rep = numpy.mean(this_rep)
        this_dam_gro = numpy.mean(this_dam_gro)
        this_dam_rep = numpy.mean(this_dam_rep)
        this_age = numpy.mean(this_age)
        this_prot = this_gro+this_rep+this_dam_gro+this_dam_rep
        this_prop_dam = (this_dam_gro+this_dam_rep)/this_prot
        this_prop_rep = this_rep/this_prot
        prot_ob.append(this_prot)
        dam_ob.append(this_prop_dam)
        rep_ob.append(this_prop_rep)
        age_ob.append(this_age)
    
    lists_a, lists_b = [prot_oa, dam_oa, rep_oa, age_oa], [prot_ob, dam_ob, rep_ob, age_ob]
    t_a1, t_b1 = t_a, t_b
    for i in range(4):
        list1, list2, list3, list4 = t_a1, lists_a[i], t_b1, lists_b[i]
        t_a, lists_a[i] = (list(x) for x in zip(*sorted(zip(list1, list2), key=itemgetter(0))))
        t_b, lists_b[i] = (list(x) for x in zip(*sorted(zip(list3, list4), key=itemgetter(0))))
    return t_a, lists_a, lists_b

main_folder = '/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/'
os.chdir(main_folder)
files = os.listdir(main_folder)
names = ['Maximum specific\ngrowth rate', 'Growth yield', 'Repair yield', 'Division radius', 'Substrate concentration', 'Substrate affinity', 'Proportional aging rate']
symbol = [r'$\mu_{max}$', r'$Y_{\mu}$', r'$Y_{r}$', r'$P_{div}$', r'$S_{bulk}$', r'$K_S$', r"$a'$"]
units = [r'h$^{-1}$', r'g g$^{-1}$', r'g g$^{-1}$', r'$\mu$m', r'g L$^{-1}$', r'g L$^{-1}$', r'h$^{-1}/\mu_G$']
values = [[1.2, 0.6, 2.4], [0.444, 0.222, 0.888], [0.8, 0.444, 1], [0.8, 0.37, 1.72], [0.003556, 0.0003556, 0.03556], [0.00234, 0.0000534, 0.108], [0.22, 0.11, 0.44]]

#mumax_1, Ymu_2, Yr_3, Pdiv_4, Sbulk_5, Ks_6, a_7
original, others = [], [[], [], [], [], [], [], []]
for f in files:
    if 'original' in f:
        if '.gif' not in f:
            original.append(f)
    else:
        for n in range(1,8):
            if str(n) in f:
                this_fol = os.listdir(main_folder+'/'+f)
                for fi in this_fol:
                    if fi != '.DS_Store' and '.gif' not in fi:
                        others[n-1].append(f+'/'+fi)

#This should be the only parameter to change to switch between growth rate and age plots
#arate=True - plotting biofilms where cells are colored by age
#arate=False - plotting biofilms where cells are colored by specific growth rate
arate=True

"""
#Start of biofilm individual plots
fig = plt.figure(figsize=(16, 12))
ax_orig = [plt.subplot(3,4,5), plt.subplot(3,4,6), plt.subplot(3,4,7)]

colbar = plt.subplot(3,16,29, frameon=False)
get_colorbar(colbar, 'autumn', 'Adaptive repair', 'left', arate=arate) #AR
colbar = plt.subplot(3,16,30, frameon=False)
get_colorbar(colbar, 'cool', 'Damage segregation', 'left', arate=arate) #ASNR
colbar = plt.subplot(3,16,31, frameon=False)
get_colorbar(colbar, 'gray', 'Solute', 'left', arate=arate) #ASNR
higher, lower = [plt.subplot(3,4,1), plt.subplot(3,4,2), plt.subplot(3,4,3)], [plt.subplot(3,4,9), plt.subplot(3,4,10), plt.subplot(3,4,11)]

for o in range(len(original)):
    plot_biofilm(ax_orig[o], main_folder+original[o], arate=arate)
    plt.sca(ax_orig[o])

for a in range(len(others)):
    ch, cl = 0, 0
    label = symbol[a]+' = '+str(values[a][0])+' '+units[a]+'\n\n'
    ax_orig[0].set_ylabel(label+r'y ($\mu$m)', fontsize=fs)
    for b in range(3):
        plt.sca(higher[b])
        plt.cla()
        plt.sca(lower[b])
        plt.cla()
    for b in range(len(others[a])):
        if 'higher' in others[a][b]:
            this_ax = higher[ch]
            plt.sca(higher[ch])
            ch += 1
            higher_label = symbol[a]+' = '+str(values[a][2])+' '+units[a]+'\n\n'
        elif 'lower' in others[a][b]:
            this_ax = lower[cl]
            plt.sca(lower[cl])
            cl += 1
            lower_label = symbol[a]+' = '+str(values[a][1])+' '+units[a]+'\n\n'
        plot_biofilm(this_ax, main_folder+others[a][b], arate=arate)
    higher[1].set_title(names[a]+'\n', fontsize=fs2+2, fontweight='bold')
    higher[0].set_ylabel(higher_label+r'y ($\mu$m)', fontsize=fs)
    lower[0].set_ylabel(lower_label+r'y ($\mu$m)', fontsize=fs)
    lower[0].set_xlabel(r'x ($\mu$m)', fontsize=fs), lower[1].set_xlabel(r'x ($\mu$m)', fontsize=fs), lower[2].set_xlabel(r'x ($\mu$m)', fontsize=fs)
    
    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    os.chdir('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/vary_parameters')
    if arate:
        plt.savefig(names[a]+'.png', dpi=600, bbox_inches='tight')
    else:
        plt.savefig(names[a]+'_growthRate.png', dpi=600, bbox_inches='tight')
plt.close()
#End of biofilm individual plots
"""
"""
#Start of log biomass ratios plots
fig = plt.figure(figsize=(16,18))

ax1, ax2, ax3, ax4 = plt.subplot2grid((42,8), (0,0), colspan=2, rowspan=9), plt.subplot2grid((42,8), (0,2), colspan=2, rowspan=9), plt.subplot2grid((42,8), (0,4), colspan=2, rowspan=9), plt.subplot2grid((42,8), (0,6), colspan=2, rowspan=9)
ax5, ax6, ax7 = plt.subplot2grid((42,8), (11,1), colspan=2, rowspan=9), plt.subplot2grid((42,8), (11,3), colspan=2, rowspan=9), plt.subplot2grid((42,8), (11,5), colspan=2, rowspan=9)

ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
labs = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
#just change this line to change between plotting mean, median, max and min
calc = 'min'

orig_t, orig_biom = [], []
for a in range(len(original)):
    #t_a, [biomass_ratio, population_ratio, growthrate_ratio] = get_biomass_ratios(main_folder, original[a])
    t_a, growthrate_ratio = get_specgrowth_ratios(main_folder, original[a], calc=calc)
    orig_t.append(t_a)
    orig_biom.append(growthrate_ratio)

for a in range(len(others)):
    ax[a].set_title(names[a])
    ax[a].set_title(labs[a], loc='left', fontweight='bold')
    s1, s2 = symbol[a]+'=', ' '+units[a]
    orig_lab, lo_lab, hi_lab = s1+str(values[a][0])+s2, s1+str(values[a][1])+s2, s1+str(values[a][2])+s2
    colors, labels = ['#D4AC0D', '#2980B9', '#B03A2E'], [lo_lab, orig_lab, hi_lab]
    ma = 0
    for b in range(len(orig_t)):
        if orig_t[b][-1] > ma: ma = orig_t[b][-1]
        ax[a].plot(orig_t[b], orig_biom[b], color=colors[1])
    ch, cl = 0, 0
    for c in range(len(others[a])):
        if 'higher' in others[a][c]:
            color = colors[2]
            ch += 1
        elif 'lower' in others[a][c]:
            color = colors[0]
            cl += 1
        #t_a, [biomass_ratio, population_ratio, growthrate_ratio] = get_biomass_ratios(main_folder, others[a][c])
        t_a, growthrate_ratio = get_specgrowth_ratios(main_folder, others[a][c], calc=calc)
        if t_a[-1] > ma: ma = t_a[-1]
        ax[a].plot(t_a, growthrate_ratio, color=color)
    handles = [Patch(facecolor=colors[c], edgecolor='k', label=labels[c]) for c in range(3)]
    handles.reverse()
    ax[a].legend(handles=handles, loc='upper left', fontsize=8)
    ax[a].plot([0,ma], [0,0], 'k--')
    if ma > 500: ma = 500
    ax[a].set_xlim([0, ma])
    ax[a].set_ylim([-2, 3])
    plt.sca(ax[a])
    if a != 0 and a != 4:
        plt.yticks([])
    ax[a].set_xlabel('Time (h)')
ax[0].set_ylabel('Log(AR/DS) specific growth rate'), ax[4].set_ylabel('Log(AR/DS) specific growth rate')

os.chdir('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/vary_parameters')
plt.savefig('Log specgrowthrate ratios '+calc+'.png', bbox_inches='tight', dpi=600)
#End of random initiation log biomass ratios

axs1a, axs1b, axs1c, axs1d, axs1e = plt.subplot2grid((42,8), (22,1), colspan=2, rowspan=4), plt.subplot2grid((42,8), (26,1), colspan=2, rowspan=4), plt.subplot2grid((42,8), (30,1), colspan=2, rowspan=4), plt.subplot2grid((42,8), (34,1), colspan=2, rowspan=4), plt.subplot2grid((42,8), (38,1), colspan=2, rowspan=4)
axs2a, axs3a = plt.subplot2grid((42,8), (22,3), colspan=2, rowspan=4), plt.subplot2grid((42,8), (22,5), colspan=2, rowspan=4)
axs2b, axs3b = plt.subplot2grid((42,8), (26,3), colspan=2, rowspan=4), plt.subplot2grid((42,8), (26,5), colspan=2, rowspan=4)
axs2c, axs3c = plt.subplot2grid((42,8), (30,3), colspan=2, rowspan=4), plt.subplot2grid((42,8), (30,5), colspan=2, rowspan=4)
axs2d, axs3d = plt.subplot2grid((42,8), (34,3), colspan=2, rowspan=4), plt.subplot2grid((42,8), (34,5), colspan=2, rowspan=4)
axs2e, axs3e = plt.subplot2grid((42,8), (38,3), colspan=2, rowspan=4), plt.subplot2grid((42,8), (38,5), colspan=2, rowspan=4)
"""
fol = os.listdir('/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/side_by_side')
fol_name = '/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/side_by_side/'
sets = ['5_Sbulk', '6_Ks', '7_aging']
"""
ax = [[axs1a, axs1b, axs1c, axs1d, axs1e], [axs2a, axs2b, axs2c, axs2d, axs2e], [axs3a, axs3b, axs3c, axs3d, axs3e]]
axs1a.set_title('Substrate concentration'), axs2a.set_title('Substrate affinity'), axs3a.set_title('Proportional aging rate')
axs1a.set_title('H', loc='left', fontweight='bold'), axs2a.set_title('I', loc='left', fontweight='bold'), axs3a.set_title('J', loc='left', fontweight='bold')

for a in range(len(sets)):
    this_ax = ax[a]
    this_ax[0].plot([0,100], [0,0], 'k--')
    for b in range(len(fol)):
        if sets[a] in fol[b] and '.gif' not in fol[b] and 'ry100' not in fol[b]:
            #t_a, [biomass_ratio, population_ratio, growthrate_ratio] = get_biomass_ratios(fol_name, fol[b])
            t_a, growthrate_ratio = get_specgrowth_ratios(fol_name, fol[b], calc=calc)
            if 'higher' in fol[b]:
                col = '#B03A2E'
            else:
                col = '#D4AC0D'
            this_ax[0].plot(t_a, growthrate_ratio, color=col)
            
            t_a, [prot_oa, dam_oa, rep_oa, age_oa], [prot_ob, dam_ob, rep_ob, age_ob] = get_prop_dam_rep(fol_name, fol[b])
            this_ax[1].plot(t_a, dam_ob, color=col) #plot adaptive repair damaged protein
            this_ax[2].plot(t_a, rep_ob, color=col) #plot adaptive repair repair protein
            
            this_ax[3].plot(t_a, dam_oa, color=col) #plot damage segregation damaged protein
            this_ax[4].plot(t_a, rep_oa, color=col) #plot damage segregation repair protein
    this_ax[0].set_xlim([0, 100]), this_ax[0].set_ylim([-2, 3])
    this_ax[1].set_xlim([0, 100]), this_ax[1].set_ylim([0, 0.5])
    this_ax[2].set_xlim([0, 100]), this_ax[2].set_ylim([0, 0.25])
    this_ax[3].set_xlim([0, 100]), this_ax[3].set_ylim([0, 1])
    this_ax[4].set_xlim([0, 100]), this_ax[4].set_ylim([-0.05, 0.05])
    if a == 0:
        this_ax[0].text(-0.2, 0.5, 'Log(AR/DS)\nspecific growth rate', rotation=90, ha='left', va='center', transform=this_ax[0].transAxes)
        this_ax[1].text(-0.28, 0.5, '\n'+r'$P_{dam}/P_{tot}$', rotation=90, ha='left', va='center', transform=this_ax[1].transAxes)
        this_ax[1].text(-0.28, 0, 'Adaptive repair', rotation=90, ha='left', va='center', transform=this_ax[1].transAxes)
        this_ax[2].text(-0.28, 0.5, '\n'+r'$P_{rep}/P_{tot}$', rotation=90, ha='left', va='center', transform=this_ax[2].transAxes)
        this_ax[3].text(-0.33, 0.5, '\n'+r'$P_{dam}/P_{tot}$', rotation=90, ha='left', va='center', transform=this_ax[3].transAxes)
        this_ax[3].text(-0.33, 0, 'Damage segregation', rotation=90, ha='left', va='center', transform=this_ax[3].transAxes)
        this_ax[4].text(-0.33, 0.5, '\n'+r'$P_{rep}/P_{tot}$', rotation=90, ha='left', va='center', transform=this_ax[4].transAxes)
    else:
        this_ax[0].set_yticks([])
        this_ax[1].set_yticks([])
        this_ax[2].set_yticks([])
        this_ax[3].set_yticks([])
        this_ax[4].set_yticks([])
    this_ax[0].set_xticks([]), this_ax[1].set_xticks([]), this_ax[2].set_xticks([]), this_ax[3].set_xticks([])
    this_ax[4].set_xlabel('Time (h)')
    s1, s2 = symbol[a+4]+'=', ' '+units[a+4]
    orig_lab, lo_lab, hi_lab = s1+str(values[a+4][0])+s2, s1+str(values[a+4][1])+s2, s1+str(values[a+4][2])+s2
    colors, labels = ['#D4AC0D', '#2980B9', '#B03A2E'], [lo_lab, orig_lab, hi_lab]
    handles = [Patch(facecolor=colors[c], edgecolor='k', label=labels[c]) for c in [0,2]]
    this_ax[0].legend(handles=handles, loc='upper left', fontsize=8)

plt.subplots_adjust(hspace=0.65)
os.chdir('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/vary_parameters')
plt.savefig('Log specgrowthrate ratios sbs '+calc+'.png', bbox_inches='tight', dpi=600)
#End of log biomass ratios of side-by-side competitions
"""
#Individual biofilm plots for side-by-side simulations
sets = ['5_Sbulk', '5_Sbulk', '6_Ks', '7_aging']
for a in range(len(sets)):
    if a != 1: continue
    fig = plt.figure(figsize=(12, 8))
    ax_high = [plt.subplot(2,3,1), plt.subplot(2,3,2)]
    ax_low = [plt.subplot(2,3,4), plt.subplot(2,3,5)]
    colbar = plt.subplot2grid((4,12), (1, 8), frameon=False, rowspan=2)
    get_colorbar(colbar, 'autumn', 'Adaptive repair', 'left', arate=arate) #AR
    colbar = plt.subplot2grid((4,12), (1, 9), frameon=False, rowspan=2)
    get_colorbar(colbar, 'cool', 'Damage segregation', 'left', arate=arate) #ASNR
    colbar = plt.subplot2grid((4,12), (1, 10), frameon=False, rowspan=2)
    get_colorbar(colbar, 'gray', 'Solute', 'left', arate=arate) #solute
    hc, lc = 0, 0
    for b in range(len(fol)):
        if a == 0:
            if 'ry100' in fol[b]: continue
        elif a == 1:
            if 'ry100' not in fol[b] and 'styrofoam' not in fol[b]: continue
        if sets[a] in fol[b] and '.gif' not in fol[b]:
            if 'higher' in fol[b] and 'styrofoam' not in fol[b]:
                ax = ax_high[hc]
                hc += 1
                higher_label = symbol[a+4]+' = '+str(values[a+4][2])+' '+units[a+4]+'\n\n'
            elif 'lower' in fol[b] or 'styrofoam' in fol[b]:
                ax = ax_low[lc]
                lc += 1
                lower_label = symbol[a+4]+' = '+str(values[a+4][1])+' '+units[a+4]+'\n\n'
            if 'styrofoam' in fol[b]:
                lower_label = 'With styrofoam\n'+symbol[a+4]+' = '+str(values[a+4][2])+' '+units[a+4]+'\n'
                higher_label = r'$Y_{r}$ = 100%'+'\n'+symbol[a+4]+' = '+str(values[a+4][2])+' '+units[a+4]+'\n'
            plot_biofilm(ax, fol_name+fol[b], arate=arate)
    ax_high[1].text(1.1, 1.05, names[a+4]+'\n', ha='center', va='bottom', fontsize=fs2+2, fontweight='bold', transform=ax_high[0].transAxes)
    ax_high[0].set_ylabel(higher_label+r'y ($\mu$m)', fontsize=fs)
    ax_low[0].set_ylabel(lower_label+r'y ($\mu$m)', fontsize=fs)
    ax_low[0].set_xlabel(r'x ($\mu$m)', fontsize=fs)
    ax_low[1].set_xlabel(r'x ($\mu$m)', fontsize=fs)
    
    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    os.chdir('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/vary_parameters')
    if a == 1: names[a+4] = names[a+4]+'_repy100'
    if arate:
        plt.savefig('Side-by-side '+names[a+4]+'.png', dpi=600, bbox_inches='tight')
    else:
        plt.savefig('Side-by-side '+names[a+4]+'_growthRate.png', dpi=600, bbox_inches='tight')
    plt.close()


calc = 'min'
ax1 = plt.subplot(111)
fol_name = '/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/side_by_side/'
ry100_high_conc = ['5_Sbulk_higher_ry100(20200701_1658)', '5_Sbulk_higher_ry100(20200701_2059)', '5_Sbulk_higher_styrofoam(20200707_1439)', '5_Sbulk_higher_styrofoam(20200707_1737)']
styrofoam = []
for a in range(4):
    t_a, [biomass_ratio, population_ratio, growthrate_ratio] = get_biomass_ratios(fol_name, ry100_high_conc[a])
    #t_a, growthrate_ratio = get_specgrowth_ratios(fol_name, ry100_high_conc[a], calc=calc)
    color = '#CF6502'
    if a > 1: color = '#019932'
    #ax1.plot(t_a, growthrate_ratio, color=color)
    ax1.plot(t_a, biomass_ratio, color=color)
    ax1.plot([0,40], [0,0], 'k--')
    ax1.set_xlim([0, 40])
    ax1.set_ylim([-2, 3])
    ax1.set_xlabel('Time (h)')
    #ax1.set_ylabel('Log(AR/DS) specific growth rate')
    ax1.set_ylabel('Log(AR/DS) biomass ratio')
handles = [Patch(facecolor='#CF6502', edgecolor='k', label=r'$Y_{r}$ = 100%'), Patch(facecolor='#019932', edgecolor='k', label='With styrofoam')]
ax1.legend(handles=handles, loc='upper left', bbox_to_anchor=(0,1))
os.chdir('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/vary_parameters')
#plt.savefig('Yr100 specific growthrate '+calc+'.png', dpi=600, bbox_inches='tight')
plt.savefig('Yr100 biomass.png', dpi=600, bbox_inches='tight')    
