#!/usr/bin/python
import numpy
import os
import sys
import toolbox_basic as basic
import toolbox_results as results
import toolbox_fitness as fitness
import toolbox_idynomics


def all_cells_growth_rate(sim_dir_path):
    biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    sim_dir_path = basic.check_path(sim_dir_path)
    attributes = {'name':'population', 'name2':'mass', 'name3':'growthRate', 'name4':'specificGrowth', 'header':'time,value,value2,value3,value4'}
    results_file_path = os.path.join(sim_dir_path, 'fitness.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    file_dir = os.path.join(sim_dir_path, 'agent_Sum')
    file_dir2 = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    basic.unzip_files(file_dir2+'.zip')
    file_list = basic.file_list(file_dir)
    file_list2 = basic.file_list(file_dir2)
    count = 0
    allMembers = [[],]*len(file_list2)
    for filename in file_list2:
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output)
        cells = species.members
        for cell in cells:
            allMembers[count] = allMembers[count] + [(cell.get_specific_growth_rate(biomass_names))]
        count += 1
    count2 = 0
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        single_result = results.SingleResult()
        species = results.SpeciesOutput(output)
        cells = species.members
        for cell in cells:
            single_result.vars['time'] = output.time
            single_result.vars['value'] = cell.vars['population']
            single_result.vars['value2'] = cell.vars['mass']
            single_result.vars['value3'] = cell.vars['growthRate']
            single_result.vars['value4'] = numpy.mean(allMembers[count2])
            result_set.members.append(single_result)
        count2 += 1
        result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    basic.rm_dir(file_dir2)
    return result_set
