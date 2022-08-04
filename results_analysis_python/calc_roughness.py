#!/usr/bin/python
from __future__ import division
import numpy
import sys
import toolbox_idynomics
import os
import toolbox_results
import toolbox_basic

def calc_roughness(sim_dir_path, iterate_number):
    sim = toolbox_idynomics.SimulationDirectory(sim_dir_path)
    agent_output = sim.get_single_iterate(iterate_number).agent_output
    nI, nJ, nK, gridres = sim.find_domain_dimensions()
    
    # Need agents coordinates and domain size
    # Matrix of agentgrid of size xdom/gridres x ydom/gridres - xdom and ydom are height and width of computational domain
    # Convert agent coordinates to agentgrid elements - 1 signifies presence
    # Define bulk and biofilm as separate bodies
    # Identify biofilm front points as agent occupied elements with at least one bulk neighbour, set equivalent elements to 1 in equivalent matrix, frontgrid    
            
    
    
    # create a grid of zeros
    agentgrid = numpy.zeros((nI, nJ), dtype=int)
    
    #set elements of the grid that have a cell in to have a value of 1 - all other grids should remain as 0
    max_x = 0
    newmax_x = 0
    for r in agent_output.get_all_cells():
        x, y = float(r.vars['locationX']), float(r.vars['locationY'])
        newx = x
        x, y = x/gridres, y/gridres
        x, y = int(x), int(y)
        if agentgrid[x, y] == 0:
            agentgrid[x, y] = 1 
        if x > max_x:
            max_x = x
        if newx > newmax_x:
            newmax_x = newx
    #set front cells to have a value of 2, making sure that they aren't on the edges where the biofilm wraps around
    for i in range(max_x+2):
        for j in range(nJ):
            if agentgrid[i, j] < 1:
                continue
            # Look at the level below, if we're not on the bottom
            if i > 0:
                if agentgrid[i-1, j] == 0:
                    agentgrid[i, j] = 2
                    continue
                # Always apply cyclic search for left/right in the level below
                if agentgrid[i-1, (j-1)%nJ] == 0:
                    agentgrid[i, j] = 2
                    continue
                if agentgrid[i-1, (j+1)%nJ] == 0:
                    agentgrid[i, j] = 2
                    continue
                # Look at the level above, if we're not at the top
            if i < nI - 1:
                if agentgrid[i+1, j] == 0:
                    agentgrid[i, j] = 2
                    continue
                # Always apply cyclic search for left/right in the level above
                if agentgrid[i+1, (j-1)%nJ] == 0:
                    agentgrid[i, j] = 2
                    continue
                if agentgrid[i+1, (j+1)%nJ] == 0:
                    agentgrid[i, j] = 2
                    continue
                # Always apply cyclic search for left/right in the same level
            if agentgrid[i, (j-1)%nJ] == 0:
                agentgrid[i, j] = 2
                continue
            if agentgrid[i,(j+1)%nJ] == 0:
                agentgrid[i,j] = 2
                continue

    frontgrid = numpy.zeros((nI), dtype=int)
    for i in range(max_x+2):
        for j in range(nJ):
            if agentgrid[i, j] == 2:
                frontgrid[i] += 1

    Cfx = numpy.zeros((nI,), dtype=float)
    Pf = 0
    num = 0
    for i in range(max_x + 1):
        temp = frontgrid[i]/nJ
        Cfx[i] = 1.0 * temp
        Pf += temp
        num += (i+1) * temp

    Xf = num/Pf
    
    Sigmaf = 0
    num2 = 0
    for j in range(nI):
        num2 += abs((j+1)-Xf)*Cfx[j]
    Sigmaf = num2/Pf
    Sigma = Sigmaf/Xf

    attributes = {'name':'Sigmaf', 'name2':'Sigma', 'name3':'Xf', 'name4':'Pf', 'name5':'height', 'header':'value,value2,value3,value4,value5'}
    results_file_path = os.path.join(sim_dir_path, 'roughness.xml')
    results_output = toolbox_results.ResultsOutput(path=results_file_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    single_result = toolbox_results.SingleResult()
    single_result.vars['value'] = Sigmaf
    single_result.vars['value2'] = Sigma
    single_result.vars['value3'] = Xf
    single_result.vars['value4'] = Pf
    single_result.vars['value5'] = max_x
    result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    #return Sigmaf, max_x

#calc_roughness(sys.argv[1])