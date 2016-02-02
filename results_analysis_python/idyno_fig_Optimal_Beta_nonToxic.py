#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting
import math
import toolbox_fitness


path = sys.argv[1]
#sim_dir1 = toolbox_idynomics.SimulationDirectory(path1)
muS = 0.3
Yr = 0.8
beta_list = []

#for iter_info in sim_dir1.get_iterate_information():
    #age = collate_mean_attributes(path1, 'age', 'agent_state', process_file_for_mean_agent_attribute,
                                                        #starting_time=240)
    #beta for toxic
    #beta = ( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1)
    #beta for non-toxic
    #beta = ( age / (1 - age))  * ( math.sqrt(Yr / muS) - 1)
    #beta_list.append(beta)




#for iter_info in sim_dir1.get_iterate_information():
    #age = float(iter_info.agent_output.species_outputs[0].calc_mean_attribute('age'))
    #beta for toxic
    #beta = ( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1)
    #beta for non-toxic
    #beta = ( age / (1 - age))  * ( math.sqrt(Yr / muS) - 1)
    #beta_list.append(beta)

aging_rates = [0.05, 0.15, 0.25]
#CE = constant environment, DE = dynamic environment/chemostat
#NT = inert/non-toxic, T = toxic
#AS = asymmetric, MS = mid-symmetric, S = symmetric
#Clegg et al. data
CE_NT_AS = [0.01, 0.02, 0.04]
CE_NT_MS = [0.01, 0.03, 0.06]
CE_NT_S = [0.01, 0.04, 0.07]
CE_T_AS = [0.03, 0.11, 0.13]
CE_T_MS = [0.03, 0.12, 0.27]
CE_T_S = [0.03, 0.13, 0.27]
DE_NT_AS = [0.03, 0.06, 0.05]
DE_NT_MS = [0.03, 0.05, 0.08]
DE_NT_S = [0.03, 0.08, 0.10]
DE_T_AS = [0.06, 0.12, 0]
DE_T_MS = [0.05, 0.14, 0.05]
DE_NT_S = [0.05, 0.14, 0.25]

calc_beta_list = CE_NT_AS

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)
axis.plot(aging_rates, calc_beta_list, 'k-')
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
fig.save(os.path.join(path, 'figures', 'optimal_beta_CE_NT_AS.png'))
