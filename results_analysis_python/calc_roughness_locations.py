#!/usr/bin/python
import numpy
import os
import sys
import toolbox_basic as basic
import toolbox_results as results
import toolbox_fitness as fitness
import toolbox_idynomics


def get_grid_size(sim_dir_path):
    sim = toolbox_idynomics.SimulationDirectory(sim_dir_path)    
    nI, nJ, nK, res = sim.find_domain_dimensions()
    grid = [nI, nJ, nK, res]
    return grid

def get_all_cells_location(sim_dir_path, starting_time=48):
    sim_dir_path = basic.check_path(sim_dir_path)
    attributes = {'name':'locations', 'header':'X,Y'}
    results_file_path = os.path.join(sim_dir_path, 'cell_locations.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    #requirements = {'family':'', 'genealogy':''}
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        if output.time >= starting_time:
            output = results.AgentOutput(path=filename)
            cells = output.get_all_cells()
            #species = results.SpeciesOutput(output)
            #cells = species.find_cells(requirements)
            for cell in cells:
                single_result = results.SingleResult()
                single_result.vars['X'] = cell.vars['locationX']
                single_result.vars['Y'] = cell.vars['locationY']
                result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    return result_set


#sim_dir_path = sys.argv[1]
#get_grid_size(sim_dir_path)
#get_all_cells_location(sim_dir_path)