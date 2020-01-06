#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import math
import numpy
import os
import toolbox_basic
import xml.etree.ElementTree as xmlTree


class Output:
    def __init__(self, path=None, root_name='idynomics'):
        if path == None or not os.path.isfile(path):
            self.path = path
            self.root = xmlTree.Element(root_name)
            self.tree = xmlTree.ElementTree(self.root)
            simulation = xmlTree.SubElement(self.root, 'simulation')
            simulation.set('iterate', '0')
            simulation.set('time', '0.0')
            simulation.set('unit', 'h')
        else:
            self.path = toolbox_basic.check_path(path)
            self.tree = toolbox_basic.get_xml_tree(self.path)
            self.root = self.tree.getroot()
        self.simulation = self.find('./simulation')
        self.iterate = int(self.simulation.attrib['iterate'])
        self.time = float(self.simulation.attrib['time'])
        self.time_unit = self.simulation.attrib['unit']
    def set_iterate(self, iterate):
        self.simulation.attrib['iterate'] = iterate
        self.iterate = iterate
    def set_time(self, time):
        self.simulation.attrib['time'] = time
        self.time = time
    def write(self, output_path=None):
        if output_path == None: output_path = self.path
        else:                   self.path = output_path
        with open(output_path, 'w') as f:
            self.tree.write(f, encoding='utf-8', xml_declaration=True)
    def find(self, search_pattern):
        out = self.tree.find(search_pattern)
        #if out == None:
        #    print('No results searching for '+search_pattern)
        return out
    def findall(self, search_pattern):
        return self.tree.findall(search_pattern)
    def display(self):
        xmlTree.dump(self.tree)


class AgentOutput(Output):
    def __init__(self, path=None, root_name='idynomics'):
        Output.__init__(self, path=path, root_name=root_name)
        if path == None:
            grid = xmlTree.SubElement(self.simulation, 'grid')
            grid.set('resolution', '0.0')
            grid.set('nI','0')
            grid.set('nJ','0')
            grid.set('nK','0')
        grid = self.find('./simulation/grid')
        self.grid_res = float(grid.attrib['resolution'])
        self.grid_nI = int(grid.attrib['nI'])
        self.grid_nJ = int(grid.attrib['nJ'])
        self.grid_nK = int(grid.attrib['nK'])
        self.three_dim = (self.grid_nK > 1)
        self.species = self.findall('./simulation/species')
        self.species_outputs = self.get_species_outputs()
    def check_single_species(self):
        if len(self.species) > 1:
            toolbox_basic.error_message('More than one species in:',self.path)
            return False
        elif len(self.species) < 1:
            toolbox_basic.error_message('No species present in:',self.path)
            return False
        else:
            return True
    def get_species_names(self):
        names = []
        for species in self.species:
            names.append(species.attrib['name'])
        return names
    def get_species_outputs(self):
        spec_list = []
        for spec_name in self.get_species_names():
            spec_list.append(SpeciesOutput(self, name=spec_name))
        return spec_list
    def get_species_by_name(self, name):
        for species in self.get_species_outputs():
            if name == species.name:
                return species
        toolbox_basic.error_message('Species %s cannot be found in'%(name), self.path)
    def get_all_cells(self):
        cell_list = []
        for spec in self.species_outputs:
            cell_list += spec.members
        return cell_list
    def add_species(self, data, header, name):
        species = xmlTree.SubElement(self.simulation, 'species')
        species.set('name', name)
        species.set('header', header)
        species.text = data
        self.species = self.findall('./simulation/species')
    def calc_total_attibute(self, attribute):
        total = 0.0
        for species in self.species_outputs:
            total += species.calc_total_attibute(attribute)
        return total


class SpeciesOutput:
    def __init__(self, agent_output, name=''):
        if not (name == ''):
            search_pattern = './simulation/species[@name="'+name+'"]'
        elif agent_output.check_single_species():
            search_pattern = './simulation/species'
        else:
            toolbox_basic.error_message('Please define which species to use in',
                                                            agent_output.path)
        self.agent_output = agent_output
        species = self.agent_output.find(search_pattern)
        self.name = species.attrib['name']
        self.attributes = species.attrib['header'].split(',')
        self.biomass_names = []
        self.members = []
        for line in species.text.translate(None,'\n').split(';'):
            if line == '': break
            variables = line.split(',')
            cell = CellOutput(self)
            for i, attribute in enumerate(self.attributes):
                cell.vars[attribute] = variables[i]
            self.members.append(cell)
    def calc_mean_specific_growth_rate(self):
        rates = self.get_specific_growth_rates()
        mean = numpy.mean(rates)
        std = numpy.std(rates)
        return mean, std
    def calc_total_specific_growth_rate(self):
        rates = self.get_specific_growth_rates()
        return sum(rates)
    def calc_total_attibute(self, attribute):
        if attribute == 'specific growth rate':
            return self.calc_total_specific_growth_rate()
        if self.attributes.count(attribute) < 1:
            toolbox_basic.error_message('Species '+self.name,
                                'does not have attribute '+attribute)
            exit()
        return sum(self.get_attribute_values(attribute))
    def calc_mean_attribute(self, attribute):
        if attribute == 'specific growth rate':
            return self.calc_mean_specific_growth_rate()
        if self.attributes.count(attribute) < 1:
            toolbox_basic.error_message('Species '+self.name,
                                'does not have attribute '+attribute)
            exit()
        values = self.get_attribute_values(attribute)
        mean = numpy.mean(values)
        std = numpy.std(values)
        return mean, std
    def change_header(self, new_header):
        search_pattern = './simulation/species[@name='+self.name+']'
        self.agent_output.species[0].attrib['header'] = new_header
        self.attributes = new_header.split(',')
    def get_header(self):
        return self.agent_output.species[0].attrib['header']
    def find_attribute_position(self, attribute):
        position = -1
        for i, x in enumerate(self.header.split(',')):
            if str(x) == str(attribute):
                position = i
                break
        if position < 0:
            msg = 'Could not find attribute "'+attribute
            msg += '" for species "'+self.name+'" in '
            toolbox_basic.error_message(msg, path)
        return position
    def find_cells(self, requirements):
        possibles = self.members
        for attribute in requirements.keys():
            requirement = str(requirements[attribute])
            possibles = [c for c in possibles if
                            (str(c.vars[attribute]) == requirement)]
        return possibles
    def get_specific_growth_rates(self):
        rates = []
        for cell in self.members:
            rates.append(cell.get_specific_growth_rate(self.biomass_names))
        return rates
    def get_attribute_values(self, attribute):
        if attribute == 'specific growth rate':
            return self.get_specific_growth_rates()
        values = []
        for cell in self.members:
            values.append(float(cell.vars[attribute]))
        return values
    def set_biomass_names(self, biomass_names):
        self.biomass_names = biomass_names
    def update_agent_output(self):
        data_script = '\n'
        for cell in self.members:
            for attribute in self.attributes[:-1]:
                data_script += str(cell.vars[attribute])+','
            data_script += str(cell.vars[self.attributes[-1]])+';\n'
        search_pattern = './simulation/species[@name='+self.name+']'
        self.agent_output.species[0].text = data_script
    def population(self):
        return len(self.members)


class CellOutput:
    def __init__(self, species):
        self.species = species.name
        self.vars = {}
        # color should be in RGB, values between 0 and 1: (r, g, b)
        self.color = None
    def get_location(self):
        x = float(self.vars['locationX'])
        y = float(self.vars['locationY'])
        if 'locationZ' in self.vars.keys():
            z = float(self.vars['locationZ'])
        else:
            z = 0.0
        return (x, y, z)
    def get_radius(self, total_radius=True):
        if total_radius:
            return float(self.vars['totalRadius'])
        else:
            return float(self.vars['radius'])
    def get_specific_growth_rate(self, biomass_names):
        growth_rate = float(self.vars['growthRate'])
        biomass = self.get_total_biomass(biomass_names)
        return growth_rate/biomass
    def get_total_biomass(self, biomass_names):
        biomass = 0.0
        for bname in biomass_names:
            biomass += float(self.vars[bname])
        return biomass
    def calc_sphere_volume(self, total_radius=True):
        #if total_radius: r = float(self.vars['totalRadius'])
        #else:         r = self.vars['radius']
        r = self.get_radius(total_radius=total_radius)
        return (4/3) * math.pi * (r**3)
    def calc_circle_area(self, total_radius=True):
        #if total_radius: r = self.vars['totalRadius']
        #else:         r = self.vars['radius']
        r = self.get_radius(total_radius=total_radius)
        return math.pi * (r**2)


class EnvOutput(Output):
    def __init__(self, path, root_name='idynomics'):
        Output.__init__(self, path, root_name=root_name)
        # If the simulation is a biofilm one, there will be a thickness element
        thickness = self.find('./simulation/thickness')
        if thickness == None:
            self.biofilm = False
        else:
            self.biofilm = True
            self.thickness_mean   = float(thickness.find('mean').text)
            self.thickness_stddev = float(thickness.find('stddev').text)
            self.thickness_max    = float(thickness.find('max').text)
        self.solutes = self.findall('./simulation/solute')
        self.solute_outputs = self.get_solute_outputs()
        
    def get_solute_outputs(self):
        sol_list = []
        for sol_name in self.get_solute_names():
            sol_list.append(SoluteOutput(self, name=sol_name))
        return sol_list
        
    def get_solute(self, solute_name):
        for solute in self.solute_outputs:
            if solute.name == solute_name:
                return solute
        toolbox_basic.error_message('Could not find solute '+solute_name,
                            'in '+self.path)

    def get_solute_names(self):
        names = []
        for solute in self.solutes:
            names.append(solute.attrib['name'])
        return names




class SoluteOutput:
    def __init__(self, env_output, name=''):
        search_pattern = './simulation/'
        if not env_output.find('./simulation/bulk') == None:
            search_pattern += 'bulk/'
        if not (name == ''):
            search_pattern += 'solute[@name="'+name+'"]'
        else:
            toolbox_basic.error_message('Please define which solute to use in',
                                                            env_output.path)
        self.env_output = env_output
        solute = env_output.find(search_pattern)
        if solute == None:
            toolbox_basic.error_message('Trouble finding solute from name:',
                                                            search_pattern)
        self.name = solute.attrib['name']
        self.unit = solute.attrib['unit']
        self.grid_res = float(solute.attrib['resolution'])
        self.grid_nI = int(solute.attrib['nI'])
        self.grid_nJ = int(solute.attrib['nJ'])
        self.grid_nK = int(solute.attrib['nK'])
        self.three_dim = (self.grid_nK > 1)
        temp = solute.text.translate(None,' ').split(';\n')
        self.values = []
        for value in temp:
            if value == '' or value == '\n': continue
            self.values.append(float(value))
    def get_concentration(self):
        return self.values[0]
    def concentration_array(self):
        self.array = numpy.array(self.values)
        if self.three_dim:
            # Older versions of iDynoMiCS included padding in the env_State
            if self.array.shape[0] == self.grid_nI*self.grid_nJ*self.grid_nK:
                new_shape = (self.grid_nI, self.grid_nJ, self.grid_nK)
                self.array = self.array.reshape(new_shape)
            else:
                new_shape = (self.grid_nI+2, self.grid_nJ+2, self.grid_nK+2)
                self.array = self.array.reshape(new_shape)
                self.array = self.array[1:-1, 1:-1, 1:-1]
        else:
            # Older versions of iDynoMiCS included padding in the env_State
            if self.array.shape[0] == self.grid_nI*self.grid_nJ:
                new_shape = (self.grid_nI, self.grid_nJ)
                self.array = self.array.reshape(new_shape)
            else:
                new_shape = (self.grid_nI+2, self.grid_nJ+2)
                self.array = self.array.reshape(new_shape)
                self.array = self.array[1:-1, 1:-1]
        return self.array


class ResultsOutput(Output):
    def __init__(self, path=None, root_name='idynomics'):
        Output.__init__(self, path=path, root_name=root_name)
        if self.get_results() == None:
            results = xmlTree.SubElement(self.simulation, 'results')
    def add_result(self, attributes, data):
        results = self.get_results()
        new_result = xmlTree.SubElement(results, 'result')
        for attribute in attributes.keys():
            new_result.set(attribute, attributes[attribute])
        new_result.text = data
        return new_result
    def get_result(self, attributes, ignore_header=False):
        attrib = attributes.copy()
        if ignore_header:
            attrib.pop('header')
        for result in self.get_results():
            result_attrib = result.attrib.copy()
            if ignore_header:
                result_attrib.pop('header')
            if toolbox_basic.are_dicts_same(result_attrib, attrib):
                return result
        return None
    def get_results(self):
        return self.find('./simulation/results')


class ResultSet:
    def __init__(self, results_output, attributes, ignore_header=False):
        self.results_output = results_output
        self.name = attributes['name']
        self.members = []
        self.result = results_output.get_result(attributes, ignore_header)
        if self.result == None:
            self.attributes = attributes['header'].split(',')
            self.result = results_output.add_result(attributes, '')
        else:
            self.attributes = self.result.attrib['header'].split(',')
            for line in self.result.text.translate(None,'\n').split(';'):
                if line == '': break
                variables = line.split(',')
                result = SingleResult()
                for i, attribute in enumerate(self.attributes):
                    result.vars[attribute] = variables[i]
                self.members.append(result)
    def update_results_output(self):
        data_script = '\n'
        for result in self.members:
            for attribute in self.attributes[:-1]:
                data_script += str(result.vars[attribute])+','
            data_script += str(result.vars[self.attributes[-1]])+';\n'
        self.result.text = data_script
        self.result.set('name', self.name)
        self.result.set('header', ','.join(self.attributes))
    def find_single_result(self, attribute_name, attribute_value):
        for result in self.members:
            if str(result.vars[attribute_name]) == str(attribute_value):
                return result
        toolbox_basic.error_message('Could not find result with '+str(attribute_name),
                            'of '+str(attribute_value))
        return None


class SingleResult:
    def __init__(self):
        self.vars = {}


'''def get_single_species(path):
    print('toolbox_results.get_single_species(path) is deprecated')
    output = AgentOutput(path)
    species = SpeciesOutput(output)
    return output, species'''
