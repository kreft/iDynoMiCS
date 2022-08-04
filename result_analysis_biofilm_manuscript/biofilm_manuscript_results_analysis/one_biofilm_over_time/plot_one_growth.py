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

folder = '/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/'
"""
#shrinking
f = 'biofilm_T_SAR_ASNR_16_cells_1(20180809_1109)'
#styrofoam
f = 'AR_DS_16cells(20180817_1932)'
"""
#proportional
f = 'biofilm_T_SAR_ASNR_16_cells_22'

species_color_dict = ['#FE9806', '#D22E03', '#d8ebf9', '#023c63']
times = os.listdir(folder+f+'/agent_State')
print(times)

for a in range(len(times)):
    if times[a][0] == 'a':
        os.chdir(folder+f+'/agent_State/')
        axis = plt.subplot(111)
        output = results.AgentOutput(path=times[a])
        time = (output.time)
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
                        cell.color=species_color_dict[0]
                    elif growthrate == 0:
                        cell.color=species_color_dict[1]
                elif species.name == 'OldieB':
                    if growthrate == 1:
                        cell.color=species_color_dict[2]
                    elif growthrate == 0:
                        cell.color=species_color_dict[3]
        toolbox_idynomics.plot_cells_2d(axis, output)
        os.chdir(folder+f+'/env_State')
        output = results.EnvOutput(path='env_State'+times[a][11:])
        for solute in output.solute_outputs:
            if solute.name == 'glucose':
                solute_output = solute
                toolbox_idynomics.solute_contour(axis, solute_output, interpolation='bicubic')
        axis.text(250, 250, 'Time: '+str(int(time))+'h', va='top', ha='right', color='#bd0303', fontsize=fs)
        axis.set_xlabel(r'$\mu$m', fontsize=fs)
        axis.set_ylabel(r'$\mu$m', fontsize=fs)
        os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/proportional/')
        plt.savefig('Growth biofilm time '+times[a][12:-5]+'.png', dpi=300, bbox_inches='tight')
        plt.close()


"""
#run this part in python3!
import os
import imageio
images = []
os.chdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/shrinking/')
filenames = os.listdir('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/shrinking/')
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('/Users/u1560915/Documents/OneDrive/Aging_of_Biofilms/Write up/paper_july_2018/one_biofilm_over_time/shrinking/movie.gif', images)
"""   
            