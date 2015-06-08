#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting


parser = OptionParser()
parser.add_option("-r", "--ResultsDir", dest="results_dir", default=os.getcwd(),
                                      help="path to results directory")
parser.add_option("-i", "--IterateNum", dest="iter_num", default=0,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-a", "--PlotAll", dest="plot_all", default=False,
                    action="store_true", help="plot all iterates, ignoring -s")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                    help="name of the solute to be plotted behind cells")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=True,
                    action="store_false", help="include a colorbar")

(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

''' TODO replace with a better way of coloring cells '''
colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)]
species_color_dict = {}
for species_name in sim.get_species_names():
    species_color_dict[species_name] = colors[0]
    colors.pop(0)



def plot(iter_info):
    nI = iter_info.agent_output.grid_nI
    nJ = iter_info.agent_output.grid_nJ
    res = iter_info.agent_output.grid_res
    # Reduce the height
    nI /= 3
    width = toolbox_plotting.mm2inch(nJ*res)
    height = toolbox_plotting.mm2inch(nI*res)
    figure = toolbox_plotting.SlideFigure(width=width, height=height)
    axis = figure.add_subplot('', 111, frameon=False)
    num_capsules = 0
    num_transconjugants = 0
    num_recipients = 0
    num_donors = 0
    for cell in iter_info.agent_output.get_all_cells():
        if cell.species == 'EPS':
            #continue
            cell.color = '0.5'
            num_capsules += 1
        elif cell.species == 'Donor':
            cell.color = 'r'
            num_donors += 1
        elif cell.species == 'Recipient':
            if int(cell.vars['copyNumber']) > 0:
                cell.color = 'g'
                num_transconjugants += 1
            else:
                cell.color = 'b'
                num_recipients += 1
        toolbox_idynomics.draw_cell_2d(axis, cell, total_radius=True)
    print('%d capsules'%(num_capsules))
    print('%d donors'%(num_donors))
    print('%d transconjugants'%(num_transconjugants))
    print('%d recipients'%(num_recipients))
    #toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
    #axis.set_title(r'Biofilm (%s g.L$^{-1}$)'%(solute_name))
    width = nJ * res
    height = nI * res
    base_thickness = 0.25*res
    axis.fill_between([0, width], [0]*2, y2=[-base_thickness]*2,
                                        facecolor='k', zorder=-1)
    axis.set_xlim(-base_thickness, width+base_thickness)
    axis.set_ylim(-base_thickness, height+base_thickness)
    lb, rt = 0.01, 0.99
    figure.subplots_adjust(left=lb, bottom=lb, right=rt, top=rt)
    figure.save(os.path.join(sim.figures_dir,
                'biofilm_%s(%d).png'%(options.solute_name, iter_info.number)))


if options.plot_all:
    min_max_concns = sim.get_min_max_concns()
    for i in sim.get_iterate_information():
        plot(i)
else:
    iter_info = sim.get_single_iterate(options.iter_num)
    plot(iter_info)
