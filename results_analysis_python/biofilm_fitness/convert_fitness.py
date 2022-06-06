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
    print 'Found path'
    """    
    attributes = {'name':'family','name1':'genealogy','name2':'generation','name3':'birthday','name4':'activeBiomassGrowth','name5':'activeBiomassRepair',
                  'name6':'inactiveBiomassGrowth','name7':'inactiveBiomassRepair','name8':'growthRate','name9':'volumeRate','name10':'locationX',
                  'name11':'locationY','name12':'locationZ','name13':'radius','name14':'totalRadius','name15':'birthX','name16':'birthY','name17':'birthZ',
                  'name18':'age','name19':'hasDied','name20':'deathX','name21':'deathY','name22':'deathZ','name23':'deathday','name24':'death', 'header':'time,value,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,value21,value22,value24'}
    """
    attributes = {'name':'population', 'name2':'mass', 'name3':'growthRate', 'name4':'specificGrowth', 'header':'time,value,value2,value3,value4'}
    for i in range(2):
        if i == 0:
            specname = 'OldieA'
        elif i == 1:
            specname = 'OldieB'
        string = str(specname)+'fitness.xml'
        results_file_path = os.path.join(sim_dir_path, string)
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
            species = results.SpeciesOutput(output, specname)
            cells = species.members
            for cell in cells:
                allMembers[count] = allMembers[count] + [(cell.get_specific_growth_rate(biomass_names))]
            count += 1
        count2 = 0
        for filename in file_list:
            output = results.AgentOutput(path=filename)
            single_result = results.SingleResult()
            species = results.SpeciesOutput(output, specname)
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
