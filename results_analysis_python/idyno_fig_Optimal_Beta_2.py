#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting
import math
import toolbox_fitness_RW

muS = 0.3
Yr = 0.8
# Set dir path
path1 = sys.argv[1]
path2 = sys.argv[2]
path3 = sys.argv[3]
sim_dir1 = toolbox_idynomics.SimulationDirectory(path1)
ages005 = toolbox_fitness_RW.calc_mean_attribute(path1, 'age', 'agent_State', toolbox_fitness_RW.process_file_for_mean_agent_attribute, starting_time=240)
ages015 = toolbox_fitness_RW.calc_mean_attribute(path2, 'age', 'agent_State', toolbox_fitness_RW.process_file_for_mean_agent_attribute, starting_time=240)
ages025 = toolbox_fitness_RW.calc_mean_attribute(path3, 'age', 'agent_State', toolbox_fitness_RW.process_file_for_mean_agent_attribute, starting_time=240)
beta = []
for age in ages005:
    #beta for toxic
    beta.append(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
for age in ages015:
    #beta for toxic
    beta.append(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
for age in ages015:
    #beta for toxic
    beta.append(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))

#beta for non-toxic
#beta.append(( age / (1 - age))  * ( math.sqrt(Yr / muS) - 1))
#beta for toxic
#beta.append(( age / (1 - age))  * ( math.sqrt((Yr / muS) * (1 / (1 - age))) - 1))
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
DE_T_AS = [0.06, 0.12]
DE_T_MS = [0.05, 0.14, 0.05]
DE_T_S = [0.05, 0.14, 0.25]

calc_beta_list = DE_T_AS
    

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)
#plot just for DE_T_AS or any others where a value takes 0
axis.plot(aging_rates[:len(calc_beta_list)], calc_beta_list, 'k-')
#axis.plot(aging_rates, calc_beta_list, 'k-')
axis.plot(0.05, beta[0], 'r*')
axis.plot(0.15, beta[2], 'r*')
axis.plot(0.25, beta[4], 'r*')
axis.set_xlim([0.05, 0.25])
axis.set_ylim([0, 1])
axis.set_xlabel('Aging rates')
axis.set_ylabel('Optimal Beta')
fig.process_subplots()
fig.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95)
#adjust for which strategy it is, and also for where the figure should be saved - could remove
#and add to e.g. path1
#path4 = sys.argv[4]
fig.save(os.path.join(path1, 'figures', 'optimal_beta_DE_T_AS.png'))
