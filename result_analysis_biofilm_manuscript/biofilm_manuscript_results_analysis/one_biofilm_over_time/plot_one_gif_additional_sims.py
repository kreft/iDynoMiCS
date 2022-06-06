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

# folder = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/SAR_ASNR/16_cells/'
# f = 'biofilm_T_SAR_ASNR_16_cells_1(20180809_1109)'
# times = os.listdir(folder+f+'/agent_State/')
# os.chdir(folder+f+'/agent_State/')

species_color_dict = {'OldieA' : 'autumn', 'OldieB' : 'cool'}

def single_time(fol, time):
    os.chdir(fol+'/agent_State/')
    axis = plt.subplot(111)
    output = results.AgentOutput(path=time)
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
    output = results.EnvOutput(path='env_State(last).xml')
    for solute in output.solute_outputs:
        if solute.name == 'glucose':
            solute_output = solute
            toolbox_idynomics.solute_contour(axis, solute_output, interpolation='bicubic')
    axis.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    axis.set_xlabel(r'$\mu$m', fontsize=fs)
    axis.set_ylabel(r'$\mu$m', fontsize=fs)
    os.mkdir(fol+'/times/')
    os.chdir(fol+'/times/')
    plt.savefig('Biofilm time '+time[12:-5]+'.png', bbox_inches='tight')
    plt.close()
    
def get_gif(folder):
    import imageio
    images = []
    biofilms = sorted(os.listdir(folder+'/times/'))
    os.chdir(folder+'/times/')
    for filename in biofilms:
        images.append(imageio.imread(filename))
    imageio.mimsave(folder+'movie_growth.gif', images, fps=3)
    
all_times = ["agent_State(0).xml","agent_State(100).xml","agent_State(200).xml","agent_State(300).xml","agent_State(400).xml","agent_State(500).xml","agent_State(600).xml","agent_State(700).xml","agent_State(800).xml","agent_State(900).xml","agent_State(1000).xml","agent_State(1100).xml","agent_State(1200).xml","agent_State(1300).xml","agent_State(1400).xml","agent_State(1500).xml","agent_State(1600).xml","agent_State(1700).xml","agent_State(1800).xml","agent_State(1900).xml","agent_State(2000).xml","agent_State(2100).xml","agent_State(2200).xml","agent_State(2300).xml","agent_State(2400).xml"]

main_folder = '/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/'
os.chdir(main_folder)
files = os.listdir(main_folder)

#mumax_1, Ymu_2, Yr_3, Pdiv_4, Sbulk_5, Ks_6, a_7
original, others = [], [[], [], [], [], [], [], []]
for f in files:
    if 'original' in f:
        original.append(f)
    else:
        for n in range(1,8):
            if str(n) in f:
                this_fol = os.listdir(main_folder+'/'+f)
                for fi in this_fol:
                    if fi != '.DS_Store':
                        others[n-1].append(f+'/'+fi)
all_folders = original
for a in range(len(others)):
    all_folders += others[a]

for f in all_folders:
    if f != 'original_random(20200603_0225)':
        continue
    for t in all_times:
        single_time(main_folder+f, t)
    get_gif(main_folder+f)
            