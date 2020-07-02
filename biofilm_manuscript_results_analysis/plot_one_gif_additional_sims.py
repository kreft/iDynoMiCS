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
import imageio
startTime = datetime.now()


fs, fs2, wspace, hspace, lb = 10, 12, 0.1, 0.4, 0.06
xticks = [0, 50, 100, 150, 200, 250]

# folder = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/SAR_ASNR/16_cells/'
# f = 'biofilm_T_SAR_ASNR_16_cells_1(20180809_1109)'
# times = os.listdir(folder+f+'/agent_State/')
# os.chdir(folder+f+'/agent_State/')

species_color_dict = {'OldieA' : 'cool', 'OldieB' : 'autumn'}

def single_time(fol, time_fn, count):
    ax = plt.subplot(111)
    os.chdir(fol+'/agent_State/')
    output = results.AgentOutput(path=time_fn)
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
    toolbox_idynomics.plot_cells_2d(ax, output)
    os.chdir(fol+'/env_State/')
    output = results.EnvOutput(path=time_fn.replace('agent', 'env'))
    for solute in output.solute_outputs:
        if solute.name == 'glucose':
            solute_output = solute
            toolbox_idynomics.solute_contour(ax, solute_output, interpolation='bicubic')
    ax.text(250, 250, 'Time: '+str(count)+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
    os.chdir(fol)
    plt.savefig(fol+'/times/'+str(count)+'.png', bbox_inches='tight', dpi=300)
    plt.close()
    return

    
def get_gif(folder):
    times_png = ["0.png","1.png","2.png","3.png","4.png","5.png","6.png","7.png","8.png","9.png","10.png","11.png","12.png","13.png","14.png","15.png","16.png","17.png","18.png","19.png","20.png","21.png","22.png","23.png","24.png",
                 "25.png","26.png","27.png","28.png","29.png","30.png","31.png","32.png","33.png","34.png","35.png","36.png","37.png","38.png","39.png","40.png","41.png","42.png","43.png","44.png","45.png","46.png","47.png","48.png"]
    #times_png = ["0.png","1.png","2.png","3.png","4.png","5.png","6.png","7.png","8.png","9.png","10.png","11.png","12.png","13.png","14.png","15.png","16.png","18.png","19.png","20.png","21.png","22.png","23.png","24.png"]
    images = []
    biofilms = times_png
    os.chdir(folder+'/times/')
    for filename in biofilms:
        try:
            images.append(imageio.imread(filename))
        except:
            time_too_high = True
    imageio.mimsave(folder+'_movie_growth.gif', images, fps=3)
    
all_times = ["agent_State(0).xml","agent_State(100).xml","agent_State(200).xml","agent_State(300).xml","agent_State(400).xml","agent_State(500).xml","agent_State(600).xml","agent_State(700).xml","agent_State(800).xml","agent_State(900).xml","agent_State(1000).xml","agent_State(1100).xml","agent_State(1200).xml","agent_State(1300).xml","agent_State(1400).xml","agent_State(1500).xml","agent_State(1600).xml","agent_State(1700).xml","agent_State(1800).xml","agent_State(1900).xml","agent_State(2000).xml","agent_State(2100).xml","agent_State(2200).xml","agent_State(2300).xml","agent_State(2400).xml",
             "agent_State(2500).xml","agent_State(2600).xml","agent_State(2700).xml","agent_State(2800).xml","agent_State(2900).xml","agent_State(3000).xml","agent_State(3100).xml","agent_State(3200).xml","agent_State(3300).xml","agent_State(3400).xml","agent_State(3500).xml","agent_State(3600).xml","agent_State(3700).xml","agent_State(3800).xml","agent_State(3900).xml","agent_State(4000).xml","agent_State(4100).xml","agent_State(4200).xml","agent_State(4300).xml","agent_State(4400).xml","agent_State(4500).xml","agent_State(4600).xml","agent_State(4700).xml","agent_State(4800).xml"]

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
for a in os.listdir('/Users/robynwright/Documents/OneDrive/Papers_writing/Aging of biofilms/Review 2/New simulations/results_files/side_by_side/'):
    all_folders.append('side_by_side/'+a)

for f in all_folders:
    if '.gif' in f:
        continue
    try:
        os.mkdir(main_folder+'/'+f+'/times/')
    except:
        do_nothing=True
    try:
        get_gif(main_folder+f)
        continue
    except:
        count = 0 
        for t in all_times:
            try:
                single_time(main_folder+f, t, count)
            except:
                time_too_high = True
            count += 1

        get_gif(main_folder+f)
            