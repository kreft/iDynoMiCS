#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

color_dict_path = os.path.join(sim_dir.figures_dir, 'color_info.txt')
if os.path.isfile(color_dict_path):
    species_color_dict = toolbox_idynomics.read_color_dict(color_dict_path)
else:
    species_color_dict = toolbox_idynomics.get_default_species_colors(sim_dir)

toolbox_idynomics.save_color_dict(species_color_dict, color_dict_path)
