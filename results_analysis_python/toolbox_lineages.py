#!/usr/bin/python
from __future__ import division
import os
import toolbox_basic as basic
import toolbox_results as results
import math
import numpy


class LineageCell(results.CellOutput):
    def __init__(self, species):
        results.CellOutput.__init__(self, species)
        self.vars['family'] = 0
        self.vars['genealogy'] = 0
        self.vars['generation'] = 0
        self.vars['initial_time'] = -1.0
        self.vars['division_time'] = -1.0
        self.vars['initial_mass'] = 1.0
        self.vars['division_mass'] = 2.0
        self.vars['growth_rate'] = 0.0
        self.sibling_mass = 0.0
    def set_growth_rate(self):
        top = math.log( self.vars['division_mass'] / self.vars['initial_mass'] )
        bottom = self.vars['division_time'] - self.vars['initial_time']
        self.vars['growth_rate'] = top / bottom


def calc_birth_generation(genealogy):
    generation = 0
    while (2**generation <= int(genealogy)):
        generation += 1
    return generation


def calc_parent_genealogy(genealogy):
    birth_generation = calc_birth_generation(genealogy)
    parent_genealogy = genealogy - 2**(birth_generation - 1)
    return parent_genealogy


def process_last_iter(sim_dir_path, max_generation, use_masses=True):
    sim_dir_path = basic.check_path(sim_dir_path)
    last_path = os.path.join(sim_dir_path, 'lastIter', 'agent_State(last).xml')
    last_path = basic.check_path(last_path)
    last_output, last_species = results.get_single_species(last_path)
    if use_masses:
        agent_State_dir = basic.unzip_files(os.path.join(sim_dir_path, 'agent_State.zip'))
        biomass_names = ['activeBiomass','inactiveBiomass']
        timestep = float(last_output.time) / float(last_output.iterate)
    # do a first run-through
    temp_cells = []
    families = []
    for cell in last_species.members:
        family = int(cell.vars['family'])
        if families.count(family) == 0:
            families.append(family)
        genealogy = int(cell.vars['genealogy'])
        birthday = float(cell.vars['birthday'])
        birth_generation = calc_birth_generation(genealogy)
        for generation in range(birth_generation, max_generation+2):
            temp = LineageCell(last_species)
            temp.vars['family'] = family
            temp.vars['genealogy'] = genealogy
            temp.vars['generation'] = generation
            temp.vars['initial_time'] = birthday
            if use_masses:
                iter = str(int( birthday / timestep ))
                iter_path = os.path.join(agent_State_dir, 
                                            'agent_State('+iter+').xml')
                iter_output, iter_species = results.get_single_species(iter_path)
                requirements = {'family':family, 'genealogy':genealogy}
                iter_cell = iter_species.find_cells(requirements)[0]
                temp.vars['initial_mass'] = iter_cell.get_total_biomass(biomass_names)
                requirements['genealogy'] = genealogy - 2**(generation-1)
                sibling = iter_species.find_cells(requirements)
                if len(sibling) == 1:
                    temp.sibling_mass = sibling[0].get_total_biomass(biomass_names)
            temp_cells.append(temp)
    print(str(len(temp_cells))+' cells in '+str(len(families))+' families')
    # then filter
    last_species.members = []
    for family in families:
        members = [c for c in temp_cells if c.vars['family'] == family]
        next_gen = [c for c in members if c.vars['generation'] == 0]
        for generation in range(max_generation + 1):
            current_gen = next_gen
            next_gen = [c for c in members if 
                                    c.vars['generation'] == generation+1]
            for cell in current_gen:
                genealogy = cell.vars['genealogy']
                # find my baby
                baby_genealogy = genealogy + 2**generation
                possibles = [c for c in next_gen if 
                                        c.vars['genealogy'] == baby_genealogy]
                if len(possibles) == 1:
                    baby = possibles[0]
                    baby_bday = baby.vars['initial_time']
                    cell.vars['division_time'] = baby_bday
                    if use_masses:
                        total_mass = baby.vars['initial_mass'] + baby.sibling_mass
                        cell.vars['division_mass'] = total_mass
                    possibles = [c for c in next_gen if 
                                            c.vars['genealogy'] == genealogy]
                    if len(possibles) == 1:
                        me_next = possibles[0]
                        me_next.vars['initial_time'] = baby_bday
                        if use_masses:
                            me_next.vars['initial_mass'] = baby.sibling_mass
                    cell.set_growth_rate()
                    last_species.members.append(cell)
    print(str(len(last_species.members))+' cells remain')
    #for cell in temp_cells.members:
    header = 'family,genealogy,generation,initial_time,division_time,'
    header += 'initial_mass,division_mass,growth_rate'
    last_species.change_header(header)
    last_species.update_agent_output()
    save_path = os.path.join(sim_dir_path, 'lineage_results.xml')
    last_output.write(output_path=save_path)


def normalise_by_generation(sim_dir_path, max_generation):
    results_path = basic.check_path(os.path.join(sim_dir_path, 'lineage_results.xml'))
    results_output = results.AgentOutput(path=results_path)
    results_species = results.SpeciesOutput(results_output)
    for generation in range(max_generation+1):
        cohort = [c for c in 
                    results_species.members if int(c.vars['generation']) == generation]
        growth_rates = [float(c.vars['growth_rate']) for c in cohort]
        mean = numpy.mean(growth_rates)
        for cell in cohort:
            cell.vars['growth_rate'] = str(float(cell.vars['growth_rate']) / mean)
    results_species.update_agent_output()
    save_path = os.path.join(sim_dir_path, 'norm_gen_lineage_results.xml')
    results_output.write(output_path=save_path)


def normalise_by_progenitor(sim_dir_path):
    results_path = basic.check_path(os.path.join(sim_dir_path, 'lineage_results.xml'))
    results_output = results.AgentOutput(path=results_path)
    results_species = results.SpeciesOutput(results_output)
    prog_growth_rates = {}
    for cell in results_species.members:
        generation = cell.vars['generation']
        if generation == '0':
            family = cell.vars['family']
            growth_rate = float(cell.vars['growth_rate'])
            prog_growth_rates[family] = growth_rate
    for cell in results_species.members:
        family = cell.vars['family']
        top = float(cell.vars['growth_rate'])
        bottom = prog_growth_rates[family]
        cell.vars['growth_rate'] = str( top / bottom )
    results_species.update_agent_output()
    save_path = os.path.join(sim_dir_path, 'norm_prog_lineage_results.xml')
    results_output.write(output_path=save_path)


def calc_averages(results_path, max_generation):
    results_path = basic.check_path(results_path)
    results_output = results.AgentOutput(path=results_path)
    results_species = results.SpeciesOutput(results_output)
    results_species.change_header(
                    'family,genealogy,generation,growth_rate,std_growth_rate')
    save_cells = []
    min_num_cells = 10000000000
    max_num_cells = 0
    for generation in range(max_generation+1):
        cohort = [c for c in results_species.members if
                        int(c.vars['generation']) == generation]
        genealogies = list(set([int(c.vars['genealogy']) for c in cohort]))
        for genealogy in genealogies:
            growth_rates = [float(c.vars['growth_rate']) for c in cohort if 
                                        int(c.vars['genealogy']) == genealogy]
            num_cells = len(growth_rates)
            min_num_cells = min(min_num_cells, num_cells)
            max_num_cells = max(max_num_cells, num_cells)
            if num_cells < 5:
                print('Only '+str(num_cells)+' in point '+str(genealogy)+', '+
                                                               str(generation))
            cell = LineageCell(results_species)
            cell.vars['genealogy'] = genealogy
            cell.vars['generation'] = generation
            cell.vars['growth_rate'] = numpy.mean(growth_rates)
            cell.vars['std_growth_rate'] = numpy.std(growth_rates)
            save_cells.append(cell)
    print('All lineages had between '+str(min_num_cells)+
                                ' and '+str(max_num_cells)+' cells')
    results_species.members = save_cells
    results_species.update_agent_output()
    dir_name = os.path.dirname(results_path)
    base_name = 'mean_'+os.path.basename(results_path)
    save_path = os.path.join(dir_name, base_name)
    results_output.write(output_path=save_path)


def plot_lineages(axis, results_path, max_generation,
                                            line_width=1.0, error_bars=True):
    results_path = basic.check_path(results_path)
    results_output = results.AgentOutput(path=results_path)
    results_species = results.SpeciesOutput(results_output)
    if error_bars and results_species.attributes.count('std_growth_rate') == 0:
        print('Could not find data for error bars')
        error_bars = False
    for cell in results_species.members:
        generation = int(cell.vars['generation'])
        if generation <= max_generation-1:
            family = int(cell.vars['family'])
            genealogy = int(cell.vars['genealogy'])
            growth_rate = float(cell.vars['growth_rate'])
            if error_bars:
                eb_x = [generation+1]
                eb_y = [growth_rate]
                eb_e = [float(cell.vars['std_growth_rate'])]
            requirements = {'family':family, 'genealogy':genealogy,
                                            'generation':generation+1}
            temp = results_species.find_cells(requirements)
            if not temp == []:
                oldie = temp[0]
                oldie_grow = float(oldie.vars['growth_rate'])
                axis.plot([generation+1, generation+2],
                          [growth_rate, oldie_grow], 
                          'r', linewidth=line_width, 
                          solid_capstyle='round', 
                          zorder=5+numpy.random.randint(5))
                if error_bars:
                    eb_x.append(generation+2)
                    eb_y.append(oldie_grow)
                    eb_e.append(float(oldie.vars['std_growth_rate']))
            requirements['genealogy'] = genealogy + 2**generation
            temp = results_species.find_cells(requirements)
            if not temp == []:
                newbie = temp[0]
                newbie_grow = float(newbie.vars['growth_rate'])
                axis.plot([generation+1, generation+2],
                          [growth_rate, newbie_grow],
                          'b', linewidth=line_width,
                          solid_capstyle='round',
                          zorder=5+numpy.random.randint(5))
                if error_bars:
                    eb_x.append(generation+2)
                    eb_y.append(newbie_grow)
                    eb_e.append(float(newbie.vars['std_growth_rate']))
            if error_bars:
                axis.errorbar(eb_x, eb_y, yerr=eb_e, ecolor='grey', fmt=None, 
                                  elinewidth=line_width/2, capsize=1, zorder=1)
    axis.set_xlabel('Generations')
    axis.set_ylabel('Normalised growth rate')
    axis.set_xlim([1, max_generation+1])


def get_first_fifteen(sim_dir_path):
    sim_dir_path = basic.check_path(sim_dir_path)
    last_path = os.path.join(sim_dir_path, 'lastIter', 'agent_State(last).xml')
    output, species = results.get_single_species(last_path)
    species.members = species.members[:15]
    for i in range(15):
        cell = species.members[i]
        cell.vars['family'] = i + 1
        cell.vars['genealogy'] = 0
        cell.vars['generation'] = 0
        cell.vars['birthday'] = 0.0
    species.update_agent_output()
    save_path = os.path.join(sim_dir_path, 'agentS.xml')
    output.write(output_path=save_path)


def setup_heterogeneous_progenitors(sim_dir_path):
    sim_dir_path = basic.check_path(sim_dir_path)
    last_path = os.path.join(sim_dir_path, 'lastIter', 'agent_State(last).xml')
    output, species = results.get_single_species(last_path)
    timestep = float(output.time) / float(output.iterate)
    agent_State_dir = basic.unzip_files(os.path.join(sim_dir_path,
                                                            'agent_State.zip'))
    save_cells = []
    for cell in species.members:
        genealogy = int(cell.vars['genealogy'])
        if genealogy == 1:
            iter = int( float(cell.vars['birthday']) / timestep )
            iter_path = os.path.join(agent_State_dir, 
                                              'agent_State('+str(iter)+').xml')
            iter_output, iter_species = results.get_single_species(iter_path)
            requirements = {'family': cell.vars['family'], 'genealogy':0}
            oldie = species.find_cells(requirements)[0]
            oldie.vars['generation'] = 0
            save_cells.append(oldie)
            requirements['genealogy'] = 1
            newbie = species.find_cells(requirements)[0]
            newbie.vars['family'] = int(newbie.vars['family']) + 15
            newbie.vars['generation'] = 0
            newbie.vars['genealogy'] = 0
            newbie.vars['birthday'] = 0.0
            save_cells.append(newbie)
    species.members = save_cells
    species.update_agent_output()
    save_path = os.path.join(sim_dir_path, 'agentS.xml')
    output.write(output_path=save_path)
