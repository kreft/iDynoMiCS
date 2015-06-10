#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

color_dict_path = os.path.join(sim_dir.figures_dir, 'color_info.txt')
if os.path.isfile(color_dict_path):
    species_color_dict = toolbox_idynomics.read_color_dict(color_dict_path)
else:
    species_color_dict = toolbox_idynomics.get_default_species_colors(sim_dir)
    toolbox_idynomics.save_color_dict(species_color_dict, color_dict_path)

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)

plot_left, plot_right, text_left = 30, 35, 38

max_time = 1
max_concn = 0

for solute_name in sim_dir.get_solute_names():
    if solute_name == 'boundaryLayer' or '-rate' in solute_name:
        continue
    time = []
    concn = []
    for iter_info in sim_dir.get_iterate_information():
        time.append(iter_info.time)
        max_time = max(max_time, iter_info.time)
        solute = iter_info.env_output.get_solute(solute_name)
        concn.append(solute.get_concentration())
        max_concn = max(max_concn, solute.get_concentration())
    axis.plot(time, concn, color=species_color_dict[solute_name])

axis.set_xlim([0, max_time])
axis.set_ylim([0, max_concn])

axis.set_xlabel('Time (h)')
axis.set_ylabel('Concentration (g/L)')

fig.process_subplots()
fig.subplots_adjust(left=0.2, bottom=0.15, right=0.98, top=0.98)
fig.save(os.path.join(path, 'figures', 'solute_concentrations.png'))
