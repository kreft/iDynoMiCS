#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

total_biomass = []
time = []
for iter_info in sim_dir.get_iterate_information():
    time.append(iter_info.time)
    total_biomass.append(iter_info.agent_output.calc_total_attibute('biomass'))
    #total_biomass.append(math.log10(iter_info.agent_output.calc_total_attibute('biomass')))

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)
axis.plot(time, total_biomass, 'k-')
#axis.plot(time, total_biomass, 'o', markerfacecolor='none', markeredgecolor='k')
#axis.plot(time, total_biomass, '+', color='0.5')
axis.set_xlim([min(time), max(time)])
#axis.set_ylim([5.5, 7.5])
axis.set_xlabel('Time (h)')
axis.set_ylabel('Total biomass (fg)')
#axis.set_ylabel(r'Log$_{10}$(Total biomass (fg))')
fig.process_subplots()
fig.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95)
fig.save(os.path.join(path, 'figures', 'total_biomass_timecourse.png'))
