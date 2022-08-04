#!/usr/bin/python
from __future__ import division
import math
import numpy
import os
from scipy import stats
import toolbox_basic as basic
import toolbox_results as results



# This is for getting the mean value of an attribute from a single simulation
# Attribute should be something like 'glucose' or 'specific growth rate'
# Output_type should something like  'agent_State' or 'env_Sum'
def calc_mean_attribute(sim_dir_path, attribute, output_type, file_process,
                                                        starting_time=2400):
    attributes = {'name':attribute, 
                    'starting_time':str(starting_time),
                    'header':'mean,std'}
    sim_dir_path = basic.check_path(sim_dir_path)
    results_file_path = os.path.join(sim_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    if len(result_set.members) > 0:
        result = result_set.members[0]
        return float(result.vars['mean']), float(result.vars['std'])
    values = []
    file_dir = os.path.join(sim_dir_path, output_type)
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    for filename in file_list:
        args = [attribute, filename]
        val = file_process(*args)
        # If output.time < starting_time, None will be returned
        if not val == None:
            values.append(val)
    single_result = results.SingleResult()
    mean_value = numpy.mean(values)
    std_value = numpy.std(values)
    single_result.vars['mean'] = mean_value
    single_result.vars['std'] = std_value
    result_set.members.append(single_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    basic.rm_dir(file_dir)
    return mean_value, std_value


# This is the file_process function for getting a mean value of an attribute 
# across the population. It should be used in conjunction with output_type 
# 'agent_State', 'agent_StateDeath', agent_Sum' or 'agent_SumDeath'
def process_file_for_mean_agent_attribute(attribute, filename,
                         biomass_names=['activeBiomassGrowth','activeBiomassRepair',
                                        'inactiveBiomassGrowth', 'inactiveBiomassRepair'], starting_time=2400):
    output = results.AgentOutput(path=filename)
    if output.time >= starting_time:
        species = results.SpeciesOutput(output)
        species.set_biomass_names(biomass_names)
        mean, std = species.calc_mean_attribute(attribute)
        return mean


# This is for collating the single simulation results of calc_mean_attribute() 
# that are together in a directory
def collate_mean_attributes(set_dir_path, attribute, output_type, file_process,
                                                        starting_time=2400):
    set_dir_path = basic.check_path(set_dir_path)
    attributes = {'name':attribute, 
                    'starting_time':str(starting_time),
                    'header':'directory,mean,std'}
    results_file_path = os.path.join(set_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    if len(result_set.members) > 0:
        return result_set
    dir_list = basic.subdir_list(set_dir_path)
    for dir_name in dir_list:
        dir_result = results.SingleResult()
        dir_result.vars['directory'] = os.path.basename(dir_name)
        mean, std = calc_mean_attribute(dir_name, attribute, output_type,
                                   file_process, starting_time=starting_time)
        dir_result.vars['mean'] = mean
        dir_result.vars['std'] = std
        result_set.members.append(dir_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    return result_set
