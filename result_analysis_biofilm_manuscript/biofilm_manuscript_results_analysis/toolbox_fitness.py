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
                         biomass_names=['activeBiomass', 'inactiveBiomass'],
                                                        starting_time=2400):
    output = results.AgentOutput(path=filename)
    if output.time >= starting_time:
        species = results.SpeciesOutput(output)
        species.set_biomass_names(biomass_names)
        mean, std = species.calc_mean_attribute(attribute)
        return mean


# This is the file_process function for getting a solute concentration.
# It should be used in conjunction with output_type 'env_State' or 'env_Sum'
def process_file_for_solute_concn(solute_name, filename,
                                                starting_time=2400):
    output = results.EnvOutput(path=filename)
    if output.time >= starting_time:
        solute = results.SoluteOutput(output, solute_name)
        return solute.get_concentration()


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


# a for asymmetric, m for mid-symmetric, s for symmetric
# Directory names should start with these
def group_dirs_by_a_m_s(set_dir_path):
    groupA = basic.subdir_list(set_dir_path, 'a*')
    groupM = basic.subdir_list(set_dir_path, 'm*')
    groupS = basic.subdir_list(set_dir_path, 's*')
    return [groupA, groupM, groupS]


def get_optima(set_dir_path, attribute, dir_group_process,
                                higher_is_fitter=True, starting_time=2400):

    results_file_path = os.path.join(set_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    attributes = {'name':'optimum '+attribute, 
                    'starting_time':str(starting_time),
                    'header':'directory,mean,std'}
    optimum_set = results.ResultSet(results_output, attributes)
    if len(optimum_set.members) > 0:
        return optimum_set
    attributes['name'] = attribute
    result_set = results.ResultSet(results_output, attributes)
    if result_set.members == []:
        basic.error_message('Collate mean attributes before assembling optima',
                                                    set_dir_path)
        exit()
    args = [set_dir_path]
    groups = dir_group_process(*args)
    for group in groups:
        single_result = results.SingleResult()
        means = {}
        for dir_name in group:
            dir_name = os.path.basename(dir_name)
            result = result_set.find_single_result('directory', dir_name)
            means[dir_name] = float(result.vars['mean'])
        if means == {}:
            continue
        if higher_is_fitter:
            optimum = max(means.values())
        else:
            optimum = min(means.values())
        fittest_dir = basic.get_key_from_dict_value(means, optimum)
        fittest = result_set.find_single_result('directory', fittest_dir)
        single_result.vars['directory'] = fittest.vars['directory']
        single_result.vars['mean'] = fittest.vars['mean']
        single_result.vars['std'] = fittest.vars['std']
        optimum_set.members.append(single_result)
    optimum_set.update_results_output()
    results_output.write(results_file_path)
    return optimum_set


def assemble_optima(superset_dir_path, attribute, dir_group_process,
                set_dir_type='*', higher_is_fitter=True, starting_time=2400):
    results_file_path = os.path.join(superset_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    attributes = {'name':'optimum '+attribute, 
                    'starting_time':str(starting_time),
                    'header':'directory,mean,std'}
    optima_set = results.ResultSet(results_output, attributes)
    if len(optima_set.members) > 0:
        return optimum_set
    set_dir_paths = basic.subdir_list(superset_dir_path, dirtype=set_dir_type)
    for set_dir_path in set_dir_paths:
        optimum_set = get_optima(set_dir_path, attribute, dir_group_process,
                higher_is_fitter=higher_is_fitter, starting_time=starting_time)
        for single_result in optimum_set.members:
            optima_set.members.append(single_result)
    optima_set.update_results_output()
    results_output.write(results_file_path)
    return optima_set


def settle_competition(set_dir_path):
    set_dir_path = basic.check_path(set_dir_path)
    attributes = {'name' : 'competition', 
                  'header' : 'speciesAwashout,speciesBwashout,numSims,pValue'}
    results_file_path = os.path.join(set_dir_path, 'results.xml')
    results_output = results.ResultsOutput(path=results_file_path)
    result_set = results.ResultSet(results_output, attributes)
    if len(result_set.members) > 0:
        return result_set
    dir_list = basic.subdir_list(set_dir_path)
    dir_result = results.SingleResult()
    washoutA = 0
    washoutB = 0
    for dir_name in dir_list:
        last_path = os.path.join(dir_name, 'lastIter', 'agent_Sum(last).xml')
        last_output = results.AgentOutput(path=last_path)
        species_names = last_output.get_species_names()
        if not len(species_names) == 2:
            basic.error_message('There should be 2 species for a competition',
                                                                     last_path)
        last_speciesA = results.SpeciesOutput(last_output,
                                                 name=species_names[0])
        populationA = float(last_speciesA.members[0].vars['population'])
        last_speciesB = results.SpeciesOutput(last_output,
                                                 name=species_names[1])
        populationB = float(last_speciesB.members[0].vars['population'])
        if populationA < 1 and populationB > 0:
            washoutA += 1
        if populationB < 1 and populationA > 0:
            washoutB += 1
        if populationA < 1 and populationB < 1:
            basic.error_message('Both species washed out in', last_path)
    dir_result.vars['speciesAwashout'] = washoutA
    dir_result.vars['speciesBwashout'] = washoutB
    dir_result.vars['numSims'] = len(dir_list)
    dir_result.vars['pValue'] = \
               stats.binom.cdf(min(washoutA,washoutB),(washoutA+washoutB),0.5)
    result_set.members.append(dir_result)
    result_set.update_results_output()
    results_output.write(results_file_path)
    return result_set


################################## Examples ##################################



def calc_mean_age(sim_dir_path):
    return calc_mean_attribute(sim_dir_path, 'age', 'agent_State',
                                         process_file_for_mean_agent_attribute)


def collate_mean_age(set_dir_path):
    return collate_mean_attributes(set_dir_path, 'age', 'agent_State',
                                     process_file_for_mean_agent_attribute)


def calc_mean_specific_growth_rate(sim_dir_path):
    return calc_mean_attribute(sim_dir_path, 'specific growth rate',
                         'agent_State',  process_file_for_mean_agent_attribute)


def collate_mean_specific_growth_rate(set_dir_path):
    return collate_mean_attributes(set_dir_path, 'specific growth rate',
                      'agent_State', process_file_for_mean_agent_attribute)


def calc_mean_glucose_concentration(sim_dir_path):
    return calc_mean_attribute(sim_dir_path, 'glucose',
                             'env_Sum',  process_file_for_solute_concn)


def collate_mean_glucose_concentration(set_dir_path):
    return collate_mean_attributes(set_dir_path, 'glucose',
                              'env_Sum', process_file_for_solute_concn)


def calc_mean_totalBiomass_concentration(sim_dir_path):
    return calc_mean_attribute(sim_dir_path, 'totalBiomass',
                             'env_State',  process_file_for_solute_concn)


def collate_mean_totalBiomass_concentration(set_dir_path):
    return collate_mean_attributes(set_dir_path, 'totalBiomass',
                              'env_State', process_file_for_solute_concn)


