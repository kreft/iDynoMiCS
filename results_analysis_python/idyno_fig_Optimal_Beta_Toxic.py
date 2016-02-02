#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting
import math


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)
muS = 0.3
Yr = 0.8

prep_ptot_list = []
beta_list = []
for iter_info in sim_dir.get_iterate_information():
    for i in 0:1000
        age = float(iter_info.agent_output.species_outputs[i].members[i].vars['age'])
        #beta for toxic
        beta = ( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1)
        #beta for non-toxic
        #beta = ( age / (1 - age))  * ( sqrt(Yr / muS) - 1)
        beta_list.append(beta)
    

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)
axis.plot(beta_list, prep_ptot_list, 'k-')
#axis.plot(beta_list, prep_ptot_list, 'r*')
#axis.plot(time, total_biomass, 'o', markerfacecolor='none', markeredgecolor='k')
#axis.plot(time, total_biomass, '+', color='0.5')
axis.set_xlim([0, 1])
axis.set_ylim([0, 1])
axis.set_xlabel('beta')
axis.set_ylabel('Prep/Ptot')
fig.process_subplots()
fig.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95)
#adjust for toxic/non-toxic
fig.save(os.path.join(path, 'figures', 'PrepPtot_optimalB_nontoxic.png'))
