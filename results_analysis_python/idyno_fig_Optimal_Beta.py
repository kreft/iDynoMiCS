#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting
import math
import toolbox_fitness_RW

# Set dir path
path = sys.argv[1]

ages = toolbox_fitness_RW.collate_mean_attributes(path, 'age', 'agent_State',
                                  toolbox_fitness_RW.process_file_for_mean_agent_attribute('age', 'agent_State'))

muS = 0.3
Yr = 0.8

beta = []
for age in ages:
    #beta for toxic
    #beta.append(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
    #beta for non-toxic
    beta.append(( age / (1 - age))  * ( sqrt(Yr / muS) - 1))

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

calc_beta_list = CE_NT_S
    

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)
axis.plot(calc_beta_list, aging_rates, 'k-')
axis.plot(beta_list , aging_rates, 'r')
axis.set_xlim([0, 1])
axis.set_ylim([0, 1])
axis.set_xlabel('Optimal Beta')
axis.set_ylabel('Aging rate')
fig.process_subplots()
fig.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95)
#adjust for toxic/non-toxic
fig.save(os.path.join(path, 'figures', 'optimal_beta_CE_NT_S.png'))
