#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting
import math


path1 = sys.argv[1]
sim_dir1 = toolbox_idynomics.SimulationDirectory(path1)
muS = 0.3
Yr = 0.8

beta_list = []
for iter_info in sim_dir1.get_iterate_information():
    age, temp = iter_info.agent_output.species_outputs[0].calc_mean_attribute('age')
    #beta for toxic
    #beta = ( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1)
    #beta for non-toxic
    beta = ( age / (1 - age))  * ( math.sqrt(Yr / muS) - 1)
    beta_list.append(beta)

with open('optimal_beta_values.txt', 'a') as f:
    f.write(beta_list)
