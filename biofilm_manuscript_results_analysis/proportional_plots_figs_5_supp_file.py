from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_results as results
import toolbox_plotting_age as toolbox_plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable
import toolbox_results
import matplotlib
import numpy
import math
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib as mpl
startTime = datetime.now()

fs, fs2, wspace, hspace, lb = 10, 12, 0.1, 0.4, 0.06
xticks = [0, 50, 100, 150, 200, 250]
    
def color_cells_age_single(axis, axis2, folder, f, species_color_dict, species_color_dict_2):
    os.chdir(folder+f+'/lastIter')
    output = results.AgentOutput(path='agent_State(last).xml')
    time = (output.time)
    for species in output.species_outputs:
        if species.members == []:
            continue
        for cell in species.members:
            name = species.name
            for (species_name, colormap_name) in species_color_dict.iteritems():
                if name == species_name:
                    norm = matplotlib.colors.Normalize(vmin=0,vmax=1)
                    colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                    age = float(cell.vars['age'])
                cell.color=m.to_rgba(age)
    toolbox_idynomics.plot_cells_2d(axis, output)
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    OldieA_activelayer = 0
    OldieB_activelayer = 0
    growthvec = []
    for species in output.species_outputs:
        if species.members == []:
            continue
        for cell in species.members:
            growthrate = cell.get_specific_growth_rate(biomass_names)
            growthrate = float(growthrate)
            growthvec.append(growthrate)
        maxspecGrowth = max(growthvec)
        interval = abs(maxspecGrowth*0.95)
    for species in output.species_outputs:
        if species.members == []:
            continue
        for cell in species.members:
            growthrate = cell.get_specific_growth_rate(biomass_names)
            growthrate = float(growthrate)
            if abs(growthrate-maxspecGrowth) <= interval:
                growthrate = 1
            else:
                growthrate = 0
            if species.name == 'OldieA': OldieA_activelayer += growthrate
            elif species.name == 'OldieB': OldieB_activelayer += growthrate
            if species.name == 'OldieA':
                if growthrate == 1: 
                    cell.color=species_color_dict_2[0]
                elif growthrate == 0:
                    cell.color=species_color_dict_2[1]
            elif species.name == 'OldieB':
                if growthrate == 1:
                    cell.color=species_color_dict_2[2]
                elif growthrate == 0:
                    cell.color=species_color_dict_2[3]
    toolbox_idynomics.plot_cells_2d(axis2, output)
    output = results.EnvOutput(path='env_State(last).xml')
    for solute in output.solute_outputs:
        if solute.name == 'glucose':
            solute_output = solute
            toolbox_idynomics.solute_contour(axis, solute_output, interpolation='bicubic')
            toolbox_idynomics.solute_contour(axis2, solute_output, interpolation='bicubic')
    axis.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    axis2.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    return
    
def get_figure():
    fig = plt.figure(figsize=(16, 10))
    ax1 = plt.subplot(241)
    ax2, ax3, ax4, ax5, ax6, ax7, ax8 = plt.subplot(242), plt.subplot(243), plt.subplot(244), plt.subplot(245), plt.subplot(246), plt.subplot(247), plt.subplot(248)
    ax1.set_title('4 cells'), ax2.set_title('8 cells'), ax3.set_title('16 cells'), ax4.set_title('32 cells')
    ax_shrink, ax_styro = [ax3, ax4, ax1, ax2], [ax7, ax8, ax5, ax6]
    for ax in ax_shrink:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    for ax in ax_styro:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    shrink_ylab, styro_ylab = 'Cells colored by age\n'+r'y ($\mu$m)', 'Cells colored by active layer\n'+r'y ($\mu$m)'
    ax1.set_ylabel(shrink_ylab)
    ax5.set_ylabel(styro_ylab)
    return fig, ax_shrink, ax_styro
    
def get_figure_colorbar(sp_col, sp_col_2, species1, species2):
    fs = 8
    print('Getting figure with colorbars')
    fig = plt.figure(figsize=(7.5, 5.625))
    ax1 = plt.subplot2grid((2,5), (0,0))
    ax2, ax3, ax4, ax5, ax6, ax7, ax8 = plt.subplot2grid((2,5), (0,1)), plt.subplot2grid((2,5), (0,2)), plt.subplot2grid((2,5), (0,3)), plt.subplot2grid((2,5), (1,0)), plt.subplot2grid((2,5), (1,1)), plt.subplot2grid((2,5), (1,2)), plt.subplot2grid((2,5), (1,3))
    axcol = [plt.subplot2grid((20, 27), (5, 22), rowspan=3, frameon=False), 
             plt.subplot2grid((20, 27), (5, 24), rowspan=3, frameon=False), 
                plt.subplot2grid((20, 27), (5, 26), rowspan=3, frameon=False),
            plt.subplot2grid((20, 27), (11, 22), rowspan=4, colspan=4)]
    for ax in axcol:
        ax.set_xticklabels(['']*10)
        ax.set_yticklabels(['']*10)
        for spine in ['right', 'left', 'top', 'bottom']:
            ax.spines[spine].set_color('none')
        ax.tick_params(top="off", bottom="off", right="off", left="off")
    fs = 6
    pad = 0.01
    side = "left"
    divider = make_axes_locatable(axcol[0])
    cax = divider.append_axes(side, size="100%", pad=pad, axisbg='none', frameon=True)
    cmap = mpl.cm.gray
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_ticks([])
    cax.text(-0.6, 0.5, 'Solute concentration', va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.1, 'No solute', va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2.5, 0.9, 'Concentrated\nsolute', va='center', ha='center', fontsize=fs-1, rotation=90)
    
    pad = 0.01
    side = "left"
    divider = make_axes_locatable(axcol[1])
    cax = divider.append_axes(side, size="100%", pad=pad, axisbg='none', frameon=True)
    cmap = sp_col[0]
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_ticks([])
    cax.text(-0.6, 0.5, 'Age '+species2[0], va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.1, 'Young cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.9, 'Old cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    
    pad = 0.01
    side = "left"
    divider = make_axes_locatable(axcol[2])
    cax = divider.append_axes(side, size="100%", pad=pad, axisbg='none', frameon=True)
    cmap = sp_col[1]
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_ticks([])
    cax.text(-0.6, 0.5, 'Age '+species2[1], va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.1, 'Young cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.9, 'Old cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    
    ax = axcol[3]
    ax.plot([0, 1], [0.0, 0.0], 'k-')
    ax.plot([0, 0], [0.95, 0.0], 'k-')
    ax.plot([1, 1], [0.95, 0.0], 'k-')
    ax.plot([0, 1], [0.95, 0.95], 'k-')
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlim([-0.02, 1.02])
    ms=3
    ax.plot(0.12, 0.8, 'o', color=sp_col_2[0], markersize=ms)
    ax.plot(0.12, 0.575, 'o', color=sp_col_2[1], markersize=ms)
    ax.plot(0.12, 0.35, 'o', color=sp_col_2[2], markersize=ms)
    ax.plot(0.12, 0.1225, 'o', color=sp_col_2[3], markersize=ms)
    ax.text(0.22, 0.8, species1[0]+'\nactive layer', va='center', fontsize=fs-1)
    ax.text(0.22, 0.575, species1[0], va='center', fontsize=fs-1)
    ax.text(0.22, 0.35, species1[1]+'\nactive layer', va='center', fontsize=fs-1)
    ax.text(0.22, 0.1225, species1[1], va='center', fontsize=fs-1)
    
    
    ax1.set_title('4 cells', fontsize=fs+2), ax2.set_title('8 cells', fontsize=fs+2), ax3.set_title('16 cells', fontsize=fs+2), ax4.set_title('32 cells', fontsize=fs+2)
    ax_age, ax_act = [ax3, ax4, ax1, ax2], [ax7, ax8, ax5, ax6]
    for ax in ax_age:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    for ax in ax_act:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    shrink_ylab, styro_ylab = 'Cells colored by age\n'+r'y ($\mu$m)', 'Cells colored by active layer\n'+r'y ($\mu$m)'
    ax1.set_ylabel(shrink_ylab)
    ax5.set_ylabel(styro_ylab)
    all_ax = [ax3, ax4, ax1, ax2, ax7, ax8, ax5, ax6]
    some_ax = [ax2, ax3, ax4, ax6, ax7, ax8]
    for ax in all_ax:
        ax.tick_params(axis='both', which='both', labelsize=fs-1, length=2)
    for ax in some_ax:
        ax.set_ylabel('')
        ax.set_yticks([])
    for ax in ax_age:
        ax.set_xlabel('')
        ax.set_xticks([])
    return fig, ax_age, ax_act, fs


def get_cols(comp_shrink):
    a_ASNR, a_SAR, a_SFR, b_ASNR, b_SAR, b_SFR = 'cool', 'autumn', 'Blues', 'cool_r', 'autumn_r', 'Blues_r'
    a_ASNR_g, a_SAR_g, a_SFR_g = ['#AED6F1', '#2874A6'], ['#FE9806', '#D22E03'], ['#D687FC', '#7D3C98']
    b_ASNR_g, b_SAR_g, b_SFR_g = ['#d8ebf9', '#023c63'], ['#fee706', '#f49605'], ['#ef87fc', '#ca2ad4']
    cols, cols_g = [], []
    species1, species2 = [], []
    for z in comp_shrink:
        sp, ind, count = ['', ''], 0, 0
        sp_g = ['', '']
        if z != 'figures':
            if z[:7] == 'control':
                count += 8
            while count < len(z):
                if z[count] == '_':
                    ind += 1
                    count += 1
                sp[ind] += z[count]
                count += 1
            new_species1, new_species2 = [], []
            for s in sp:
                if s == 'ASNR':
                    new_species1.append('Damage\nsegregation')
                    new_species2.append('Damage segregation')
                elif s == 'SAR':
                    new_species1.append('Adaptive\nrepair')
                    new_species2.append('Adaptive repair')
                elif s == 'SFR':
                    new_species1.append('Fixed\nrepair')
                    new_species2.append('Fixed repair')
            species1.append(new_species1)
            species2.append(new_species2)
            if sp[0] == sp[1]:
                if sp[0] == 'ASNR':
                    cols.append([a_ASNR, b_ASNR])
                    cols_g.append([a_ASNR_g, b_ASNR_g])
                elif sp[0] == 'SAR':
                    cols.append([a_SAR, b_SAR])
                    cols_g.append([a_SAR_g, b_SAR_g])
                elif sp[0] == 'SFR':
                    cols.append([a_SFR, b_SFR])
                    cols_g.append([a_SFR_g, b_SFR_g])
            else:
                for x in range(len(sp)):
                    if sp[x] == 'ASNR':
                        sp[x] = a_ASNR
                        sp_g[x] = a_ASNR_g
                    if sp[x] == 'SAR':
                        sp[x] = a_SAR
                        sp_g[x] = a_SAR_g
                    if sp[x] == 'SFR':
                        sp[x] = a_SFR
                        sp_g[x] = a_SFR_g
                cols.append(sp)
                cols_g.append(sp_g)
    return cols, cols_g, species1, species2
    
folder_prop = '/Volumes/Robyn_W_2/july_2018/paper/biofilms_proportional/'
comp_prop = sorted(os.listdir(folder_prop))
cols, cols_g, species1, species2 = get_cols(comp_prop)
comp_prop.sort()

for a in range(len(comp_prop)):
    species_color_dict = {'OldieA' : cols[a][0], 'OldieB' : cols[a][1]}
    species_color_dict_2 = [cols_g[a][0][0], cols_g[a][0][1], cols_g[a][1][0], cols_g[a][1][1]]
    folder_comp_prop = folder_prop+comp_prop[a]+'/'
    cells_prop = sorted(os.listdir(folder_comp_prop))
    files = [[], [], [], []]
    folders = []
    for b in range(len(cells_prop)):
        this_folder = folder_comp_prop+cells_prop[b]+'/'
        folders.append(this_folder)
        for c in os.listdir(this_folder):
            if c[-9:] != '.DS_Store':
                files[b].append(c)
            files[b].sort()
    reps = len(files[0])
    #for b in range(1):
    for b in range(reps):
        if a == 0 or a == 1:
            continue
        elif a == 2 and b < 38:
            continue
        fig, ax_age, ax_gro, fs = get_figure_colorbar(cols[a], species_color_dict_2, species1[a], species2[a])
        for c in range(len(files)):
            cells_prop_folder = folders[c]
            sim = files[c][b]
            print cells_prop_folder+sim
            color_cells_age_single(ax_age[c], ax_gro[c], cells_prop_folder, sim, species_color_dict, species_color_dict_2)
            print(datetime.now() - startTime)
        print('Saving figure')
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/proportional_biofilms/')
        fig.subplots_adjust(wspace=0.15, hspace=-0.6)
        plt.savefig(comp_prop[a]+'_'+str(b+1)+'.png', bbox_inches='tight', dpi=600)
        plt.close()
        print(datetime.now() - startTime)
print(datetime.now() - startTime)