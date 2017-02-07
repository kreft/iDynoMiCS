#!/usr/bin/python
import numpy
import os
import sys
import toolbox_basic as basic
import toolbox_results as results
import toolbox_fitness_RW as fitness

#This is the part that is used for Fig4A
def build_life_history(sim_dir_path, attribute1, attribute2, attribute3, attribute4, attribute5, attribute6, cell_id=(1,0),
                        biomass_names=['activeBiomassGrowth', 'activeBiomassRepair',
                                       'inactiveBiomassGrowth', 'inactiveBiomassRepair']):
    sim_dir_path = basic.check_path(sim_dir_path)
    #this must be value and not value1 for the first value so that the plotting
    #script can also plot older results at the same time
    attributes = {'name':attribute1, 'name2':attribute2, 'name3':attribute3, 'name4':attribute4, 'name5':attribute5, 'name6':attribute6, 'header':'time,value,value2,value3,value4,value5,value6'}
    results_file_path = os.path.join(sim_dir_path, 'results_for_Prep_Ptot.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    if len(result_set.members) > 0:
        return result_set.members
    requirements = {'family':str(cell_id[0]), 'genealogy':str(cell_id[1])}
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output)
        cells = species.find_cells(requirements)
        single_result = results.SingleResult()
        single_result.vars['time'] = output.time
        if len(cells) == 0:
            continue
        cell = cells[0]
        single_result.vars['value'] = cell.vars[attribute1]
        single_result.vars['value2'] = cell.vars[attribute2]
        single_result.vars['value3'] = cell.vars[attribute3]
        single_result.vars['value4'] = cell.vars[attribute4]
        single_result.vars['value5'] = cell.vars[attribute5]
        single_result.vars['value6'] = cell.vars[attribute6]
        result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    return result_set