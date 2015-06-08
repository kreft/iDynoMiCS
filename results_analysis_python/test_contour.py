#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import numpy
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


parser = OptionParser()
parser.add_option("-r", "--ResultsDir", dest="results_dir",
                  default=os.getcwd(), help="path to results directory")
parser.add_option("-i", "--IterateNum", dest="iter_num", default=-1,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                    help="name of the solute to be plotted behind cells")

(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

iter_info = sim.get_single_iterate(options.iter_num)

save_name = 'biofilm_'+options.solute_name

num_digits = len(str(sim.get_last_iterate_number()))

colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)]
species_color_dict = {}
script = 'Species name\t\tRGB color\n'
for species_name in sim.get_species_names():
    script += species_name+'\t\t'+str(colors[0])+'\n'
    species_color_dict[species_name] = colors[0]
    colors.pop(0)
toolbox_idynomics.color_cells_by_species(iter_info.agent_output,
                                                            species_color_dict)


figure = toolbox_plotting.SlideFigure()
axis = figure.add_subplot('', 111, frameon=False, projection='3d')
axis.view_init(azim=45, elev=35)

solute_output = toolbox_results.SoluteOutput(iter_info.env_output,
                                                     name=options.solute_name)

res = solute_output.grid_res
width = solute_output.grid_nJ * res
height = solute_output.grid_nI * res
depth = solute_output.grid_nK * res
extent = [0, width, 0, height]
array = solute_output.concentration_array()

min_concn = -1.0
max_concn = numpy.max(array)
'''
cs = axis.imshow(array,
            interpolation='nearest', origin='lower', cmap='gray',
            extent=extent, zorder=0)
'''

levels = numpy.linspace(0.0, max_concn, 128)

arrayX = array[:, :, 0]
Y, Z = numpy.meshgrid(numpy.linspace(0, depth, solute_output.grid_nK),
                          numpy.linspace(0, height, solute_output.grid_nI))
cs = axis.contourf(arrayX, Y, Z, zdir='x', cmap='gray', offset=0,
                               zorder=-1, levels=levels)#, vmin=min_concn, vmax=max_concn)
arrayY = array[:, 0, :]
X, Z = numpy.meshgrid(numpy.linspace(0, width, solute_output.grid_nJ),
                          numpy.linspace(0, height, solute_output.grid_nI))
cs = axis.contourf(X, arrayY, Z, zdir='y', cmap='gray', offset=0,#depth,
                               zorder=-1, levels=levels)#, vmin=min_concn, vmax=max_concn)
#arrayZ = array[1][:][:]
arrayZ = numpy.zeros([17, 17])
X, Y = numpy.meshgrid(numpy.linspace(0, width, solute_output.grid_nJ),
                          numpy.linspace(0, depth, solute_output.grid_nK))
cs = axis.contourf(X, Y, arrayZ, zdir='z', cmap='gray', offset=0,
                               zorder=-1, levels=levels)#, vmin=min_concn, vmax=max_concn)

toolbox_idynomics.plot_cells_3d(axis, iter_info.agent_output)

axis.set_xlim([0, width])
axis.set_ylim([0, depth])
axis.set_zlim([0, height])

#axis.set_xlabel('x')
#axis.set_ylabel('y')
#axis.set_zlabel('z')

#toolbox_plotting.make_colorbar(axis, cs, side="left")
#figure.fig.colorbar(cs)

figure.save('blah.png')