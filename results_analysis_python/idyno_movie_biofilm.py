#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_basic
import toolbox_idynomics


parser = OptionParser()
parser.add_option("-f", "--FrameRate", dest="frame_rate", default=24,
                        type="int", help="number of images per second")
parser.add_option("-r", "--ResultsDir", dest="results_dir",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                        help="name of the solute to be plotted behind cells")
(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

save_name = 'biofilm_'+options.solute_name

file_list = toolbox_basic.file_list(sim.figures_dir, save_name+'*.png')
num_digits = len(str(len(file_list)))

temp_dir = os.path.join(os.path.abspath(sim.movies_dir), 'temp')
toolbox_basic.make_dir(temp_dir)
for i, filename in enumerate(file_list):
    save_num = str(i)
    link_name = 'img' + (num_digits-len(save_num))*'0' + save_num + '.png'
    link_name = os.path.join(temp_dir, link_name)
    toolbox_basic.make_symbolic_link(filename, link_name)


cmd = "ffmpeg -framerate "+str(options.frame_rate)+"  -i '"
cmd += os.path.join(temp_dir, "img%"+str(num_digits)+"d.png'")
cmd += " -pix_fmt yuv420p -r 24  '"
cmd += os.path.join(os.path.abspath(sim.movies_dir), save_name+".mp4'")
os.system(cmd)

toolbox_basic.rm_dir(temp_dir)
