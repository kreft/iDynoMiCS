#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 12:42:45 2020

@author: robynwright
"""

import matplotlib.pyplot as plt
import os
import toolbox_basic as basic
import toolbox_idynomics
import toolbox_results as results
import matplotlib

concs = []
mumax = 1.2
conc = []
Kss = [0.0000534, 0.00234, 0.108]
colors = ['#D4AC0D', '#2980B9', '#B03A2E']
labels = [r'K$_S$='+str(a)+r'g L$^{-1}$' for a in Kss]
ma = 1#0.03556

a = 0
while a < ma:
    conc.append(a)
    a += 0.0000001


plt.figure(figsize=(6,5))
ax1a = plt.subplot2grid((8,2), (0,0), rowspan=3)
ax1b = plt.subplot2grid((8,2), (0,1), rowspan=3)

for a in range(len(Kss)):
    Ks = Kss[a]
    this_plot = []
    for S in conc:
        muS = (mumax*S)/(Ks+S)
        this_plot.append(muS)
    ax1a.plot(conc, this_plot, color=colors[a], label=labels[a])
    ax1b.plot(conc, this_plot, color=colors[a], label=labels[a])
    #plt.semilogx()
ax1b.legend(loc='upper left', bbox_to_anchor=(1, 1.04))
ax1a.set_ylabel(r'$\mu$(S)')
#ax1.set_xlabel(r'$S_{bulk}$ (g L$^{-1}$)')
ax1a.set_xlim([0, 0.01])
ax1b.set_xlim([0, 1])
plt.sca(ax1b)
plt.yticks([0, 0.5, 1], ['', '', ''])
plt.xticks([0, 0.5, 1], ['0', '0.5', '1'])
plt.sca(ax1a)
plt.xticks([0, 0.005, 0.01], ['0', '0.005', '0.01'])

ax2 = plt.subplot2grid((8,1), (3,0), rowspan=5)

fi = '/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/side_by_side/5_Sbulk_higher(20200625_0842)'

os.chdir(fi+'/agent_State')
output = results.AgentOutput(path='agent_State(1000).xml')
time = (output.time)
species_color_dict = {'OldieA' : 'cool', 'OldieB' : 'autumn'}
xlims = {'OldieA':[25, 75], 'OldieB':[175, 225]}
biomass_names=['activeBiomassGrowth', 'activeBiomassRepair','inactiveBiomassGrowth', 'inactiveBiomassRepair']
colors, growth_rates, xloc, yloc = [], [], [], []
for species in output.species_outputs:
    if species.members == []:
        continue
    for cell in species.members:
        name = species.name
        x, y = float(cell.vars['locationX']), float(cell.vars['locationY'])
        """
        if x > 120:
            lims = xlims[name]
            if y >= lims[0] and y <= lims[1]:
                adding_cell = True
            else:
                continue
        else:
            continue
        """
        lims = xlims[name]
        for (species_name, colormap_name) in species_color_dict.items():
            if name == species_name:
                colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                norm = matplotlib.colors.Normalize(vmin=0,vmax=1)
                m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                age = float(cell.vars['age'])
                cell.color=m.to_rgba(age)
        spec_growth = float(cell.get_specific_growth_rate(biomass_names))
        colors.append(cell.color)
        growth_rates.append(spec_growth)
        xloc.append(x)
        yloc.append(y)
os.chdir(fi+'/env_State')
output = results.EnvOutput(path='env_State(1000).xml')
for solute in output.solute_outputs:
    if solute.name != 'glucose': continue
    array = solute.concentration_array()
solute = list(array)
concs = []
for a in range(len(colors)):
    x, y = int(xloc[a]/4), int(yloc[a]/4)
    conc = solute[x][y]
    concs.append(conc)
    ax2.scatter(conc, growth_rates[a], color=colors[a], s=20, edgecolors='gray', linewidths=0.1)
ax2.set_xlim([0, 0.006])
ax2.set_ylabel(r'$\mu$(S)')
ax2.set_xlabel(r'$S$ (g L$^{-1}$)')

ax1a.set_title('A', loc='left', fontweight='bold', fontsize=12)
ax1b.set_title('B', loc='left', fontweight='bold', fontsize=12)
ax2.text(.01,.9,'C', horizontalalignment='left', transform=ax2.transAxes, fontweight='bold', fontsize=12)

plt.subplots_adjust(hspace=1, wspace=0.1)
plt.savefig('/Users/robynwright/Documents/GitHub/iDynoMiCS_1.5/biofilm_manuscript_results_analysis/additional_plots/spec_growth.png', dpi=600, bbox_inches='tight')