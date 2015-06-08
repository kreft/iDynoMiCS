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
y, y_diff = 100, -60

max_time = 1
max_abundance = 1

for species_name in sim_dir.get_species_names():
    time = []
    abundance = []
    for iter_info in sim_dir.get_iterate_information():
        time.append(iter_info.time)
        max_time = max(max_time, iter_info.time)
        species = iter_info.agent_output.get_species_by_name(species_name)
        abundance.append(species.population())
        max_abundance = max(max_abundance, species.population())
    #print species_name
    #print time
    #print abundance
    axis.plot(time, abundance, color=species_color_dict[species_name])
    #axis.plot([plot_left, plot_right], [y]*2, color=species_color_dict[species_name])
    #axis.text(text_left, y, species_name, va='center', ha='left')
    #y += y_diff

axis.set_xlim([0, max_time])
axis.set_ylim([0, max_abundance])

axis.set_xlabel('Time (h)')
axis.set_ylabel('Number of cells')

fig.process_subplots()
fig.subplots_adjust(left=0.2, bottom=0.15, right=0.98, top=0.98)
fig.save(os.path.join(path, 'figures', 'total_species_abundance.png'))
