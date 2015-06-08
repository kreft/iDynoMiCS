#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


parser = OptionParser()
parser.add_option("-a", "--AllIter", dest="all_iter", default=False,
                  action="store_true", help="plot all iterates, ignoring -i")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=False,
                            action="store_true", help="include a colorbar")
parser.add_option("-d", "--DotsPerInch", dest="dpi", default=300, type="int",
                                              help="output figure resolution")
parser.add_option("-e", "--FileExt", dest="file_ext", default=".png",
                                      help="file extension for figure output")
parser.add_option("-f", "--FrameOn", dest="frameon", default=False,
                        action="store_true", help="turn the figure frame on")
parser.add_option("-F", "--FigureType", dest="figure_type", default=None,
                         help="type of figure to use. Default is 'Slide'")
parser.add_option("-H", "--Height", dest="height", default=0,
                    type="int", help="figure height in inches")
parser.add_option("-i", "--IterateNum", dest="iter_num", default=0,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-I", "--IMax", dest="i_max", default=0,
                        type="int", help="maximum height to plot")
parser.add_option("-r", "--ResultsDir", dest="results_dir",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                        help="name of the solute to be plotted behind cells")
parser.add_option("-S", "--Substratum", dest="substratum", default=True,
                  action="store_false", help="turn off drawing the substratum")
parser.add_option("-t", "--TimeOn", dest="timeon", default=False,
                        action="store_true", help="record the time in figures")
parser.add_option("-T", "--TitleOn", dest="titleon", default=False,
                        action="store_true", help="turn the figure title on")
parser.add_option("-W", "--Width", dest="width", default=0,
                    type="int", help="figure width in inches")
parser.add_option("-z", "--ZeroColorBar", dest="zero_color", default=False,
                    action="store_true",
                    help="forces the lower limit of the color bar to zero")
(options, args) = parser.parse_args()


sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

if not options.file_ext[0] == '.':
    options.file_ext = '.'+options.file_ext

save_name = 'biofilm_'+options.solute_name

num_digits = len(str(sim.get_last_iterate_number()))
    
color_dict_path = os.path.join(sim.figures_dir, 'color_info.txt')
if os.path.isfile(color_dict_path):
    species_color_dict = toolbox_idynomics.read_color_dict(color_dict_path)
else:
    species_color_dict = toolbox_idynomics.get_default_species_colors(sim)
    toolbox_idynomics.save_color_dict(species_color_dict, color_dict_path)

nI, nJ, nK, res = sim.find_domain_dimensions()
if options.i_max > 0:
    nI = options.i_max

counter = 0

if options.figure_type == None:
    if options.height > 0: height = options.height
    else: height = toolbox_plotting.mm2inch(nI * res)
    if options.width > 0: width = options.width
    else: width = toolbox_plotting.mm2inch(nJ * res)
    figure = toolbox_plotting.SlideFigure(width=width, height=height)
else:
    script = "figure = toolbox_plotting."+options.figure_type+"Figure("
    if nI > 2*nJ:
        script += "height='double'"
    elif nJ > 2*nI:
        script += "double_column=True, height='single'"
    else:
        script += "double_column=True, height='double'"
    script += ")"
    try:
        exec(script)
    except:
        print 'Could not make figure!'
        print script

def plot(iter_info, min_max_concns):
    axis = figure.add_subplot('', 111, frameon=options.frameon)
    toolbox_idynomics.color_cells_by_species(iter_info.agent_output,
                                                           species_color_dict)
    toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
    if options.substratum:
        axis.fill_between([0, nJ*res], [0]*2, y2=[-res]*2, color='k', zorder=-1)
    figure.subplots_adjust(left=0.01, bottom=0.01, right=0.9, top=0.9)
    figure.inset_axes()
    if not options.solute_name == "none":
        #print min_max_concns[options.solute_name]
        solute_output = toolbox_results.SoluteOutput(iter_info.env_output,
                                                     name=options.solute_name)
        cs = toolbox_idynomics.solute_contour(axis, solute_output,
                            concn_range=min_max_concns[options.solute_name],
                                                    interpolation='bicubic')
        if options.color_bar:
            toolbox_plotting.make_colorbar(axis, cs)
    if options.titleon:
        axis.set_title(r'Biofilm (%s g L$^{-1}$)'%(options.solute_name))
    if options.frameon:
        axis.set_xlabel('x')
        axis.set_ylabel('y')
    if options.timeon:
        axis.text(0.0, nI*res, 'Time: %d h'%(int(iter_info.time)),
                  va='top', ha='left')
    axis.set_xlim(0, nJ * res)
    axis.set_ylim(-res, nI * res)
    save_num = str(iter_info.number)
    #save_num = str(counter)
    #counter += 1
    save_num = (num_digits - len(save_num))*'0' + save_num
    save_path = os.path.join(sim.figures_dir,
                                     save_name+'_'+save_num+options.file_ext)
    figure.save(save_path, dpi=options.dpi)


if options.all_iter:
    min_max_concns = sim.get_min_max_concns()
    if options.zero_color:
        min_max_concns[options.solute_name][0] = 0.0
    for i in sim.get_iterate_numbers():
        iter_info = sim.get_single_iterate(i)
        plot(iter_info, min_max_concns)
else:
    iter_info = sim.get_single_iterate(options.iter_num)
    min_max_concns = iter_info.get_min_max_concns()
    if options.zero_color:
        min_max_concns[options.solute_name][0] = 0.0
    plot(iter_info, min_max_concns)
