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

"""
fs, fs2, wspace, hspace, lb = 10, 12, 0.1, 0.4, 0.06
xticks = [0, 50, 100, 150, 200, 250]

folder = '/Volumes/Robyn_W_2/july_2018/paper/shrinking_biofilms/SAR_ASNR/16_cells/'
f = 'biofilm_T_SAR_ASNR_16_cells_1(20180809_1109)'

species_color_dict = {'OldieA' : 'autumn', 'OldieB' : 'cool'}
times = os.listdir(folder+f+'/agent_State/')
os.chdir(folder+f+'/agent_State/')

for a in range(len(times)):
    os.chdir(folder+f+'/agent_State/')
    axis = plt.subplot(111)
    output = results.AgentOutput(path=times[a])
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
    os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/shrinking/')
    plt.savefig('Biofilm time '+times[a][12:-5]+'.png', bbox_inches='tight')
    plt.close()
"""
biofilms_shrinking = ['Biofilm time 240.png','Biofilm time 480.png','Biofilm time 720.png','Biofilm time 960.png','Biofilm time 1200.png','Biofilm time 1440.png','Biofilm time 1680.png','Biofilm time 1920.png','Biofilm time 2160.png','Biofilm time 2400.png','Biofilm time 2640.png','Biofilm time 2880.png','Biofilm time 3120.png','Biofilm time 3360.png','Biofilm time 3600.png','Biofilm time 3638.png']
folder_shrinking = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/shrinking/'
biofilms_styrofoam = ['Biofilm time 240.png','Biofilm time 480.png','Biofilm time 720.png','Biofilm time 960.png','Biofilm time 1200.png','Biofilm time 1323.png','Biofilm time last.png']
folder_styrofoam = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/styrofoam/'
biofilms_proportional = ['Biofilm time 200.png','Biofilm time 580.png','Biofilm time 980.png','Biofilm time 1240.png','Biofilm time 2180.png','Biofilm time 2600.png','Biofilm time 3380.png','Biofilm time 3400.png','Biofilm time 3800.png','Biofilm time 4280.png','Biofilm time 4500.png','Biofilm time 4900.png','Biofilm time 5080.png','Biofilm time 5700.png']
biofilms_proportional_growth = ['Growth biofilm time 200.png','Growth biofilm time 580.png','Growth biofilm time 980.png','Growth biofilm time 1240.png','Growth biofilm time 2180.png','Growth biofilm time 2600.png','Growth biofilm time 3380.png','Growth biofilm time 3400.png','Growth biofilm time 3800.png','Growth biofilm time 4280.png','Growth biofilm time 4500.png','Growth biofilm time 4900.png','Growth biofilm time 5080.png','Growth biofilm time 5700.png']
folder_proportional = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/proportional/'
#run this part in python3!
biofilms = biofilms_proportional_growth
folder = folder_proportional
import imageio
images = []
os.chdir(folder)
for filename in biofilms:
    images.append(imageio.imread(filename))
imageio.mimsave(folder+'movie_growth.gif', images, fps=3)
 
            