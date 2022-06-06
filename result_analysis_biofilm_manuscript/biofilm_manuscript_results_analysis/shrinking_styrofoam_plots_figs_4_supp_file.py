#shrinking and styrofoam biofilms
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
    
def color_cells_age_single(axis, folder, f, species_color_dict):
    os.chdir(folder+f+'/lastIter')
    output = results.AgentOutput(path='agent_State(last).xml')
    time = (output.time)
    for species in output.species_outputs:
        if species.members == []:
            continue
        for cell in species.members:
            name = species.name
            for (species_name, colormap_name) in species_color_dict.items():
                if name == species_name:
                    norm = matplotlib.colors.Normalize(vmin=0,vmax=1)
                    colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                    age = float(cell.vars['age'])
                cell.color=m.to_rgba(age)
    toolbox_idynomics.plot_cells_2d(axis, output)
    output = results.EnvOutput(path='env_State(last).xml')
    for solute in output.solute_outputs:
        if solute.name == 'glucose':
            solute_output = solute
            toolbox_idynomics.solute_contour(axis, solute_output, interpolation='bicubic')
    axis.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    return
    
def get_figure_colorbar(sp_col, species2):
    fs = 8
    print('Getting figure with colorbars')
    fig = plt.figure(figsize=(7.5, 5.625))
    ax1 = plt.subplot2grid((2,5), (0,0))
    ax2, ax3, ax4, ax5, ax6, ax7, ax8 = plt.subplot2grid((2,5), (0,1)), plt.subplot2grid((2,5), (0,2)), plt.subplot2grid((2,5), (0,3)), plt.subplot2grid((2,5), (1,0)), plt.subplot2grid((2,5), (1,1)), plt.subplot2grid((2,5), (1,2)), plt.subplot2grid((2,5), (1,3))
    axcol = [plt.subplot2grid((20, 27), (8, 22), rowspan=4, frameon=False), 
             plt.subplot2grid((20, 27), (8, 24), rowspan=4, frameon=False), 
                plt.subplot2grid((20, 27), (8, 26), rowspan=4, frameon=False)]
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
    cax = divider.append_axes(side, size="100%", pad=pad, facecolor='none', frameon=True)
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
    cax = divider.append_axes(side, size="100%", pad=pad, facecolor='none', frameon=True)
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
    cax = divider.append_axes(side, size="100%", pad=pad, facecolor='none', frameon=True)
    cmap = sp_col[1]
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_ticks([])
    cax.text(-0.6, 0.5, 'Age '+species2[1], va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.1, 'Young cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    cax.text(2, 0.9, 'Old cells', va='center', ha='center', fontsize=fs-1, rotation=90)
    
    
    ax1.set_title('4 cells', fontsize=fs+2), ax2.set_title('8 cells', fontsize=fs+2), ax3.set_title('16 cells', fontsize=fs+2), ax4.set_title('32 cells', fontsize=fs+2)
    ax_age, ax_act = [ax3, ax4, ax1, ax2], [ax7, ax8, ax5, ax6]
    for ax in ax_age:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    for ax in ax_act:
        ax.set_xlabel(r'x ($\mu$m)', fontsize=fs)
        ax.set_ylabel(r'y ($\mu$m)', fontsize=fs)
    age_ylab, act_ylab = 'Cells without styrofoam\n'+r'y ($\mu$m)', 'Cells with styrofoam\n'+r'y ($\mu$m)'
    ax1.set_ylabel(age_ylab)
    ax5.set_ylabel(act_ylab)
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
    a_ASNR, a_SAR, a_SFR, b_ASNR, b_SAR, b_SFR = plt.cm.cool, plt.cm.autumn, plt.cm.Blues, plt.cm.cool_r, plt.cm.autumn_r, plt.cm.Blues_r
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
                    new_species1.append('damage\nsegregation')
                    new_species2.append('damage segregation')
                elif s == 'SAR':
                    new_species1.append('adaptive\nrepair')
                    new_species2.append('adaptive repair')
                elif s == 'SFR':
                    new_species1.append('fixed\nrepair')
                    new_species2.append('fixed repair')
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

folder_shrink = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/'
folder_styro = '/Volumes/Robyn_W_2/july_2018/paper/styrofoam_competitions_new/'
comp_shrink, comp_styro = sorted(os.listdir(folder_shrink)), sorted(os.listdir(folder_styro))
cols, cols_g, species1, species2 = get_cols(comp_shrink)
species_color_dict = {'OldieA' : 'autumn', 'OldieB' : 'cool'}
comp_shrink.sort()
comp_styro.sort()


count = 0
for a in range(1):
#for a in range(len(comp_shrink)):
    if comp_shrink[a] != 'figures':
        species_color_dict = {'OldieA' : cols[count][0], 'OldieB' : cols[count][1]}
        folder_comp_shrink = folder_shrink+comp_shrink[a]+'/'
        folder_comp_styro = folder_styro+comp_styro[a]+'/'
        cells_shrink, cells_styro = sorted(os.listdir(folder_comp_shrink)), sorted(os.listdir(folder_comp_styro))
        reps = len(os.listdir(folder_comp_shrink+cells_shrink[0]))
        #for b in range(reps):
        for b in range(5):
            #b += 1
            if b > 4:
                continue
            if b == 0:
                continue
            fig, ax_shrink, ax_styro, fs = get_figure_colorbar(cols[a], species2[a])
            for c in range(len(cells_shrink)):
            #for c in range(1):
                species_color_dict = {'OldieA' : cols[count][0], 'OldieB' : cols[count][1]}
                cell_shrink_folder = folder_comp_shrink+cells_shrink[c]+'/'
                sims = sorted(os.listdir(cell_shrink_folder))
                print(cell_shrink_folder+sims[b])
                color_cells_age_single(ax_shrink[c], cell_shrink_folder, sims[b], species_color_dict)
                if comp_shrink[a] == 'ASNR_SFR':
                    species_color_dict = {'OldieA' : cols[count][1], 'OldieB' : cols[count][0]}
                cell_styro_folder = folder_comp_styro+cells_shrink[c]+'/'
                sims = sorted(os.listdir(cell_styro_folder))
                print(cell_styro_folder+sims[b])
                color_cells_age_single(ax_styro[c], cell_styro_folder, sims[b], species_color_dict)
                print(datetime.now() - startTime)
            print(datetime.now() - startTime)
            os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/shrinking_styro_biofilms/')
            fig.subplots_adjust(wspace=0.15, hspace=-0.6)
            plt.savefig(comp_shrink[a]+'_'+str(b+1)+'.png', bbox_inches='tight', dpi=600)
            plt.close()
            print(datetime.now() - startTime)
        count += 1


"""
f = 'biofilm_T_SAR_SFR_32_cells_1(20180810_0245)'
folder = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/SAR_SFR/32_cells/'
os.chdir(folder)
species_color_dict = {'OldieA' : 'autumn', 'OldieB' : 'cool'}
fig = plt.figure(figsize=(5, 5))
ax1 = plt.subplot(111)
color_cells_age_single(ax1, folder, f, species_color_dict)
ax1.set_xlabel(r'x ($\mu$m)', fontsize=fs)
os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/shrinking_styro_biofilms/')
#fig.savefig(f+' age.png', bbox_inches='tight', dpi=600)
#plt.close()
"""
print(datetime.now() - startTime)
