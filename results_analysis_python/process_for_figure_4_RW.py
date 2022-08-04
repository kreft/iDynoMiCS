#!/usr/bin/python
import numpy
import os
import sys
import toolbox_basic as basic
import toolbox_results as results
import toolbox_fitness as fitness


def build_population_structure(sim_dir_path, attribute, 
                        bins=numpy.linspace(0, 0.6, 121), starting_time=2400,
                        biomass_names=['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'activeBiomassRepair']):
    attributes = {'name':attribute+' population structure',
                    'starting_time':str(starting_time),
                    'header':'bin,frequency'}
    sim_dir_path = basic.check_path(sim_dir_path)
    results_file_path = os.path.join(sim_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    if len(result_set.members) > 0:
        return result_set.members
    attr_values = []
    total_pop = 0.0
    file_dir = os.path.join(sim_dir_path, 'agent_State')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        if output.time >= starting_time:
            species = results.SpeciesOutput(output)
            species.set_biomass_names(biomass_names)
            attr_values.extend(species.get_attribute_values(attribute))
            total_pop += species.population()
    hist, bin_edges = numpy.histogram(attr_values, bins)
    for i in range(len(hist)):
        single_result = results.SingleResult()
        single_result.vars['bin'] = str(bins[i+1])+'-'+str(bins[i])
        single_result.vars['frequency'] = str(float(hist[i])/total_pop)
        result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    return result_set


# cell_id is (family, genealogy)
def build_life_history(sim_dir_path, attribute, cell_id=(1,0),
                        biomass_names=['activeBiomass', 'inactiveBiomass']):
    sim_dir_path = basic.check_path(sim_dir_path)
    attributes = {'name':attribute, 'header':'time,value'}
    results_file_path = os.path.join(sim_dir_path, 'results.xml')
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
        cell = cells[0]
        single_result = results.SingleResult()
        single_result.vars['time'] = output.time
        if attribute == 'specific growth rate':
            single_result.vars['value'] = \
                        cell.get_specific_growth_rate(biomass_names)
        else:
            single_result.vars['value'] = cell.vars[attribute]
        result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    return result_set


