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
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import matplotlib as mpl

parser = OptionParser()
'''
parser.add_option("-a", "--AllIter", dest="all_iter", default=False,
                  action="store_true", help="plot all iterates, ignoring -i")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=False,
                            action="store_true", help="include a colorbar")
                            '''
parser.add_option("-d", "--DotsPerInch", dest="dpi", default=300, type="int",
                                              help="output figure resolution")
parser.add_option("-e", "--FileExt", dest="file_ext", default=".png",
                                      help="file extension for figure output")
parser.add_option("-f", "--FrameOn", dest="frameon", default=True,
                        action="store_true", help="turn the figure frame on")
'''
parser.add_option("-F", "--FigureType", dest="figure_type", default=None,
                         help="type of figure to use. Default is 'Slide'")
parser.add_option("-H", "--Height", dest="height", default=0,
                    type="int", help="figure height in inches")
                    '''
parser.add_option("-i", "--IterateNum", dest="iter_num", default=1441,
                    type="int", help="number of the iterate to be plotted")
'''
parser.add_option("-I", "--IMax", dest="i_max", default=0,
                        type="int", help="maximum height to plot")
                        '''
parser.add_option("--rr", "--ResultsDirr", dest="results_dirr",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("--Rr", "--ResultsDiRr", dest="results_diRr",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("--RR", "--ResultsDiRR", dest="results_diRR",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="glucose",
                        help="name of the solute to be plotted behind cells")
parser.add_option("-S", "--Substratum", dest="substratum", default=True,
                  action="store_false", help="turn off drawing the substratum")
'''
parser.add_option("-t", "--TimeOn", dest="timeon", default=False,
                        action="store_true", help="record the time in figures")
parser.add_option("-T", "--TitleOn", dest="titleon", default=False,
                        action="store_true", help="turn the figure title on")
'''
parser.add_option("--TitleScript", help="script of the title",
                   dest="titlescript", default=None)
'''
parser.add_option("-W", "--Width", dest="width", default=0,
                    type="int", help="figure width in inches")
parser.add_option("-z", "--ZeroColorBar", dest="zero_color", default=False,
                    action="store_true",
                    help="forces the lower limit of the color bar to zero")
'''
parser.add_option("--sp", dest="speciesnumber", default=1, type="int", 
                  help="how many species in the biofilm, script allows for 2")
parser.add_option("--N", dest="numplots", default=1, type="int",
                  help="allows more than one biofilm to be plotted in the same figure")
(options, args) = parser.parse_args()

gs = gridspec.GridSpec(3,1)
def get_axes(numplots, plot):
    if numplots == 1:
        ax0 = plt.subplot2grid((4,1), (0, 0), colspan=3)
        ax1 = plt.subplot2grid((4,1), (0, 3), colspan=1, sharey=ax0)
    elif numplots == 2:
        if plot == 0:
            ax0 = plt.subplot2grid((8,1), (0, 0), colspan=3)
            ax1 = plt.subplot2grid((8,1), (0, 3), colspan=1, sharey=ax0)
        elif plot == 1:
            ax0 = plt.subplot2grid((8,1), (0, 4), colspan=3)
            ax1 = plt.subplot2grid((8,1), (0, 7), colspan=1, sharey=ax0)
    elif numplots == 3:
        if plot == 0:
            ax0 = plt.subplot2grid((12,1), (0, 0), colspan=3)
            ax1 = plt.subplot2grid((12,1), (0, 3), colspan=1, sharey=ax0)
        elif plot == 1:
            ax0 = plt.subplot2grid((12,1), (0, 4), colspan=3)
            ax1 = plt.subplot2grid((12,1), (0, 7), colspan=1, sharey=ax0)
        elif plot == 2:
            ax0 = plt.subplot2grid((12,1), (0, 8), colspan=3)
            ax1 = plt.subplot2grid((12,1), (0, 11), colspan=1, sharey=ax0)
    else:
        print ("Don't know which axis to use")
    axis = [ax0, ax1]
    return axis

numplots = options.numplots
if numplots == 1:
    sim1 = toolbox_idynomics.SimulationDirectory(options.results_dirr)
if numplots == 2:
    sim1 = toolbox_idynomics.SimulationDirectory(options.results_dirr)
    sim2 = toolbox_idynomics.SimulationDirectory(options.results_diRr)
if numplots == 3:
    sim1 = toolbox_idynomics.SimulationDirectory(options.results_dirr)
    sim2 = toolbox_idynomics.SimulationDirectory(options.results_diRr)
    sim3 = toolbox_idynomics.SimulationDirectory(options.results_diRR)

if not options.file_ext[0] == '.':
    options.file_ext = '.'+options.file_ext

save_name = 'biofilm_'+str(options.numplots)
num_digits = len(str(sim1.get_last_iterate_number()))

nI, nJ, nK, res = sim1.find_domain_dimensions()


species_color_dict = {'OldieA' : 'cool', 'OldieB' : 'Greens', 'OldieC' : 'summer'}
def color_cells_by_age_and_species(agent_output, species_color_dict):   
    for species in agent_output.species_outputs:
        if species.members == []:
            continue
        print('Colouring %d %s cells %s'%(len(species.members),
                            species.name, species_color_dict[species.name]))
        for cell in species.members:
            name = species.name
            for (species_name, colormap_name) in species_color_dict.iteritems():
                if name == species_name:
                    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
                    colormap = matplotlib.cm.get_cmap(colormap_name, 256)
                    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
                    age = cell.vars['age']
                    age = float(age)
                cell.color = m.to_rgba(age)
                
def colorbar(axis, position="bottom"):
    cmap = matplotlib.cm.get_cmap('cool', 256)
    #a = numpy.array([[0,1]])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cbar = mpl.colorbar.ColorbarBase(axis, cmap=cmap, norm=norm, orientation='horizontal')
    
    return cbar

agegrid = [0]*nJ
locationX = [0]*nJ
def get_agegrid():
    print('Getting cell ages')
    for r in iter_info.agent_output.get_all_cells():
        x = float(r.vars['locationX'])
        x = int(x)
        locationX[x] = x
        age = float(r.vars['age'])
        if age == 0:
            continue
        else:
            agegrid[x] += age
            agegrid[x] /= 2
    maximum = max(locationX)+1
    agegrid2 = [0]*maximum
    locationX2 = [0]*maximum
    for a in range(maximum):
        agegrid2[a] = agegrid[a]
        locationX2[a] = locationX[a]
    return agegrid2, locationX2


def plot(iter_info, min_max_concns, plot):    
    axis = get_axes(numplots, plot)
    biofilm = axis[0]
    ageplot = axis[1]
    color_cells_by_age_and_species(iter_info.agent_output, species_color_dict)
    toolbox_idynomics.plot_cells_2d(biofilm, iter_info.agent_output)
    ages, locationX = get_agegrid()
    
    #colorbar
    if numplots == 1:
        cbar = colorbar(biofilm)
    elif numplots == 2:
        cbar = colorbar(biofilm)
    elif numplots == 3 and plot == 1:
        cbar = colorbar(biofilm)
    
    ageplot.plot(ages, locationX)
    ageplot.set_yticks([])
    ageplot.set_xlim(-0.1,1)
    ageplot.set_xticks([0, 0.5, 1.0])
    ageplot.set_xlabel('Mean Age')
    
    
    biofilm.fill_between([0, nJ*res], [0]*2, y2=[-res]*2, color='k', zorder=-1)
    lb = 0.01

    #if options.frameon:
    lb = 0.06
    #process_lines()

    if not options.solute_name == "none":
        solute_output = toolbox_results.SoluteOutput(iter_info.env_output,
                                                     name=options.solute_name)
        cs = toolbox_idynomics.solute_contour(biofilm, solute_output,
                            concn_range=min_max_concns[options.solute_name],
                                                    interpolation='bicubic')
    if not options.titlescript is None:
        axis.set_title(options.titlescript)
    biofilm.set_xlabel(r'x ($\mu$m)')
    biofilm.set_ylabel(r'y ($\mu$m)')
    if plot > 0:
        biofilm.set_yticks([])
    biofilm.set_xlim(0, nJ * res)
    biofilm.set_ylim(-res, nI * res)    
    
'''   
   if plot == numplots-1:
        figure.subplots_adjust(left=lb, bottom=lb, right=0.9, top=0.9)
        figure.inset_axes()
        save_num = str(iter_info.number)
        #save_num = str(counter)
        #counter += 1
        save_num = (num_digits - len(save_num))*'0' + save_num
        save_path = os.path.join(sim.figures_dir, save_name+'_'+save_num+options.file_ext)
        figure.save(save_path, dpi=options.dpi)
 '''

for j in range(numplots):
    if j == 0:
        sim = sim1
    elif j == 1:
        sim = sim2
    elif j == 2:
        sim = sim3
    else:
        sim = sim1
    iter_info = sim.get_single_iterate(options.iter_num)
    min_max_concns = iter_info.get_min_max_concns()
    plot(iter_info, min_max_concns, j)


plt.show()