#/usr/bin/python
from __future__ import division
from __future__ import with_statement
import math
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy
from numpy import mean as amean
import os
import re
from scipy.spatial import Delaunay
from scipy.spatial import KDTree
from scipy.stats.mstats import gmean as gmean
from scipy.stats.mstats import hmean as hmean
import toolbox_basic
import toolbox_plotting
import toolbox_results_new as toolbox_results
import toolbox_schematic
import xml.etree.ElementTree as xmlTree


class Simulation:
    def __init__(self, root_dir, find_protocol=True):
        self.root_dir = toolbox_basic.check_path(root_dir)
        if not find_protocol: return
        # Set up the space from the protocol file.
        protocol_path = toolbox_basic.find_protocol_file_path(self.root_dir)
        self.protocol_tree = toolbox_basic.get_xml_tree(protocol_path)
        space_params = self.protocol_tree.findall('./space/param')
        for param in space_params:
            name, text = param.attrib['name'], param.text
            if name == 'wrapping': self.wrapping = (text == 'true')
            elif name == 'length': self.side_length = int(float(text))
            elif name == 'nDims':  self.is_3d = (int(float(text)) == 3)
        marks = self.protocol_tree.findall('./mark')
        for mark in marks:
            for param in mark:
                if param.attrib['name'] == 'value':
                    value = int(param.text)
                if param.attrib['name'] == 'number':
                    number = int(param.text)
            if value == 2:
                self.consumers = number
            else:
                self.producers = number
        rs = self.protocol_tree.find("./process/param[@name='randomSeed']")
        self.random_seed = int(rs.text)
    def copy_mesh_files(self, detail_level, sif_name):
        self.get_detail_level_mesh_dir(detail_level)
        run_dir = self.get_run_dir(detail_level, sif_name)
        run_mesh_dir = os.path.join(run_dir, 'mesh')
        toolbox_basic.copy_dir(self.detail_level_mesh_dir, run_mesh_dir)
    def get_concn_array(self, detail_level, sif_name):
        self.get_run_dir(detail_level, sif_name)
        grid_length = detail_level * self.side_length
        array_path = os.path.join(self.run_dir, 'concn_array')
        if os.path.isfile(array_path):
            concn_array = toolbox_basic.load_array(array_path)
        else:
            result_file_path = os.path.join(self.run_dir, 'case.result')
            result_file_path = toolbox_basic.check_path(result_file_path)
            # This isn't quite correct! Without wrapping, Elmer skips some nodes
            num_nodes = (grid_length + 1)**2
            array_shape = (grid_length + 1,)*2
            with open(result_file_path, 'Ur') as f:
                last_lines = f.readlines()[-num_nodes:]
            concn_array = numpy.array([float(line) for line in last_lines])
            concn_array = numpy.reshape(concn_array, array_shape)
            toolbox_basic.save_array(concn_array, array_path)
        return concn_array
    def get_consume_produce_functions(self, detail_level, sif_name):
        self.get_sif_file_path(detail_level, sif_name)
        with open(self.sif_file_path) as f:
            lines = f.readlines()
        regex = re.compile('\$ function consume.*')
        cons_line = [line for line in lines if re.match(regex, line)][0]
        min_index = cons_line.index('min(0.0')
        close_index = cons_line.index(') }')
        consume_rate = 'def consume_rate(c):\n\treturn min(0.0, %s)' \
                                        %(cons_line[min_index+8:close_index])
        regex = re.compile('\$ function produce.*')
        prod_line = [line for line in lines if re.match(regex, line)][0]
        min_index = prod_line.index('max(0.0')
        close_index = prod_line.index(') }')
        produce_rate = 'def produce_rate(c):\n\treturn max(0.0, %s)' \
                                        %(prod_line[min_index+8:close_index])
        exec consume_rate
        exec produce_rate
        return consume_rate, produce_rate
    def get_detail_level_dir(self, detail_level):
        name = 'detail_level_%d'%(detail_level)
        self.detail_level_dir = os.path.join(self.root_dir, name)
        toolbox_basic.make_dir(self.detail_level_dir)
        return self.detail_level_dir
    def get_detail_level_mesh_dir(self, detail_level):
        self.get_detail_level_dir(detail_level)
        self.detail_level_mesh_dir = os.path.join(self.detail_level_dir, 'mesh')
        toolbox_basic.make_dir(self.detail_level_mesh_dir)
        return self.detail_level_mesh_dir
    def get_detail_level_results(self, detail_level, read_only=False):
        self.get_detail_level_dir(detail_level)
        cells_path = os.path.join(self.detail_level_dir, 'cell_locations.xml')
        cells_file = SimulationResultsFile(path=cells_path, read_only=read_only)
        return cells_file
    def get_combined_results(self, detail_level, sif_name):
        dl_results = self.get_detail_level_results(detail_level, read_only=True)
        run_results = self.get_run_results(detail_level, sif_name, read_only=True)
        for e in run_results.events:
            for o in dl_results.events:
                if o.position() == e.position():
                    for attrib in o.vars.keys():
                        e.vars[attrib] = o.vars[attrib]
        return run_results
    def get_rate_array(self, detail_level, sif_name):
        self.get_run_dir(detail_level, sif_name)
        array_path = os.path.join(self.run_dir, 'rate_array')
        if not os.path.isfile(array_path):
            self.get_run_results(detail_level, sif_name)
        return toolbox_basic.load_array(array_path)
    def get_run_results(self, detail_level, sif_name, read_only=False):
        self.get_run_dir(detail_level, sif_name)
        rates_path = os.path.join(self.run_dir, 'cell_rates.xml')
        rates_file = SimulationResultsFile(path=rates_path, read_only=read_only)
        if rates_file.events == []:
            grid_length = self.side_length * detail_level
            rates_file = self.get_detail_level_results(detail_level)
            rates_file.path = rates_path
            rates_file.setup_ranges(grid_length, self.wrapping)
            rate_array = numpy.zeros((grid_length,)*2, dtype=numpy.float)
            concn_array = self.get_concn_array(detail_level, sif_name)
            consume_rate, produce_rate = \
                    self.get_consume_produce_functions(detail_level, sif_name)
            rates_file.calc_rates_from_concn_array(concn_array, consume_rate,
                                            produce_rate, rate_array=rate_array)
            array_path = os.path.join(self.run_dir, 'rate_array')
            toolbox_basic.save_array(rate_array, array_path)
            head = 'mark,x,y'
            if self.is_3d: head += ',z'
            head += ',rate,amean_surf_concn'
            rates_file.set_event_list_header(head)
            rates_file.set_concn_rate_results(concn_array, rate_array)
            rates_file.write()
        return rates_file
    def get_run_dir(self, detail_level, sif_name):
        self.get_detail_level_dir(detail_level)
        self.run_dir = os.path.join(self.detail_level_dir, sif_name)
        toolbox_basic.make_dir(self.run_dir)
        return self.run_dir
    def get_sif_file_path(self, detail_level, sif_name):
        self.get_run_dir(detail_level, sif_name)
        self.sif_file_path = os.path.join(self.run_dir, sif_name+'.sif')
        return self.sif_file_path
    def make_start_file(self, detail_level, sif_name):
        self.get_run_dir(detail_level, sif_name)
        file_path = os.path.join(self.run_dir, 'ELMERSOLVER_STARTINFO')
        with open(file_path, 'w') as f:
            f.write('''%s\n1''' %(sif_name+'.sif'))
    def make_mesh_files(self, biomass_array, detail_level):
        self.get_detail_level_mesh_dir(detail_level)
        grid_length = self.side_length * detail_level
        num_elements = (grid_length)**2
        num_nodes = (grid_length+1)**2
        ### Make mesh.header
        header_path = os.path.join(self.detail_level_mesh_dir, 'mesh.header')
        # The 2 denotes dimensions, 202's are boundaries, 404's are elements.
        with open(header_path, 'w') as f:
            f.write('%d\t%d\t%d\t\n2\t\n202\t%d\t\n404\t%d\t\n\t'
                    %(num_nodes, num_elements, num_elements, num_elements, num_elements))
        ### Make mesh.nodes
        text = ''
        for i in range(num_nodes):
            # Shouldn't this take account of detail_level?
            (y, x) = divmod(i, (grid_length+1))
            # Consider changing this line to
            #text += '%d -1 %.1f %.1f 0.0\n' %(i+1, x, y)
            text += str(i+1)+' -1 '+str(x)+' '+str(y)+' 0.0\n'
        nodes_path = os.path.join(self.detail_level_mesh_dir, 'mesh.nodes')
        with open(nodes_path, 'w') as f:
            f.write(text)
        ### Make mesh.elements
        text = ''
        counter = 0
        for (i, j), body in numpy.ndenumerate(biomass_array):
            counter += 1
            n1 = (j+1) + (i*(grid_length+1))
            n2 = n1 + 1
            n3 = n2 + (grid_length+1)
            n4 = n3 - 1
            text += '%d %d 404 %d %d %d %d \n' %(counter, body, n1, n2, n3, n4)
        elements_path = os.path.join(self.detail_level_mesh_dir, 'mesh.elements')
        with open(elements_path, 'w') as f:
            f.write(text)
        ### Make mesh.boundary
        text = ''
        counter = 0
        # Along the bottom of the array (x=max) from left (y=0) to right (y=max).
        e_base = grid_length*(grid_length - 1) + 1
        n_base = grid_length*(grid_length + 1) + 1
        for i in range(grid_length):
            counter += 1
            element = e_base + i
            node = n_base + i
            text += '%d 1 %d 0 202 %d %d \n' %(counter, element, node, node+1)
        # Down the left of the array (y=0), from top (x=0) to bottom (x=max).
        n_base = grid_length + 1
        for i in range(grid_length):
            counter += 1
            element = (i*grid_length) + 1
            node = 1 + i*n_base
            text += '%d 2 %d 0 202 %d %d \n' %(counter, element, node, node+n_base)
        # Along the top of the array (x=0) from left (y=0) to right (y=max).
        for i in range(grid_length):
            counter += 1
            text += '%d 3 %d 0 202 %d %d \n' %(counter, i+1, i+1, i+2)
        # Down the left of the array (y=max), from top (x=0) to bottom (x=max).
        n_base = grid_length + 1
        for i in range(grid_length):
            counter += 1
            element = (i+1)*grid_length
            node = (i+1)*n_base
            text += '%d 4 %d 0 202 %d %d \n' %(counter, element, node+n_base, node)
        boundary_path = os.path.join(self.detail_level_mesh_dir, 'mesh.boundary')
        with open(boundary_path, 'w') as f:
                f.write(text)
    def set_up_population(self, detail_level):
        grid_length = self.side_length * detail_level
        cells_file = self.get_detail_level_results(detail_level)
        # If cells_file.events is empty then the detail level directory has
        # probably only just been created, and this file does not yet exist.
        if cells_file.events == []:
            cells_file.set_space_parameters(wrapping=self.wrapping,
                                is_3d=self.is_3d, side_length=self.side_length)
            bio_array_path = os.path.join(self.detail_level_dir, 'bio_array')
            self.bio_array = numpy.ones((grid_length,)*2, dtype=numpy.int)
            last_path = os.path.join(self.root_dir, 'lastIter',
                                                     'event_location_last.xml')
            last_file = SimulationResultsFile(path=last_path, read_only=True)
            cells_file.copy_event_list(last_file)
            cells_file.set_up_population(detail_level, grid_length,
                                                self.wrapping, self.bio_array)
            toolbox_basic.save_array(self.bio_array, bio_array_path)
            self.make_mesh_files(self.bio_array, detail_level)
            # Finally, update and save the 'cell_locations.xml' file.
            head = 'mark,x,i_min,i_max,y,j_min,j_max'
            if self.is_3d: head += ',z,k_min,k_max'
            cells_file.set_event_list_header(head)
            cells_file.write()
        else:
            cells_file.setup_ranges(grid_length, self.wrapping)
        return cells_file
    def calc_amean_surf_concn(self, concn_array):
        for event in self.event_list:
            event.calc_amean_surf_concn(concn_array)
    def plot_concn_array(self, axis, detail_level, sif_name, set_as_white=None, plot_cs=True):
        array = self.get_concn_array(detail_level, sif_name)
        extent = [-0.5/detail_level, self.side_length + 0.5/detail_level]*2
        bottom_red, top_red     = 0.1, 0.7
        bottom_green, top_green = 0.6, 0.0
        bottom_blue, top_blue   = 0.1, 0.5
        mid_point = 0.5
        if not set_as_white == None:
            max_val, min_val = numpy.max(array), numpy.min(array)
            up_diff, down_diff = max_val - set_as_white, set_as_white - min_val
            max_diff, total_diff = max(up_diff, down_diff), max_val - min_val
            up_rel_diff, down_rel_diff = up_diff/max_diff, down_diff/max_diff
            mid_point = down_diff/total_diff
        cdict = {'red':   ((0, bottom_red, bottom_red),
                           (mid_point, 1, 1),
                           (1, top_red, top_red)),
                 'green': ((0, bottom_green, bottom_green),
                           (mid_point, 1, 1),
                           (1, top_green, top_green)),
                 'blue':  ((0, bottom_blue, bottom_blue),
                           (mid_point, 1, 1),
                           (1, top_blue, top_blue))}
        my_cmap = \
               matplotlib.colors.LinearSegmentedColormap('my_cmap', cdict, 255)
        cs = axis.imshow(array, interpolation='nearest', origin='lower',
                                                   extent=extent, cmap=my_cmap)
        axis.set_xlim(0.0, self.side_length), axis.set_ylim(0.0, self.side_length)
        if plot_cs:
            cbar = toolbox_plotting.make_colorbar(axis, cs, fontsize=8)
            return cbar
        else:
            return cs
    def plot_rate_array(self, axis, detail_level, sif_name):
        array = self.get_rate_array(detail_level, sif_name)
        extent = [0.0, self.side_length]*2
        max_val = numpy.max(abs(array))
        cdict = {'red':   ((0, 0, 0), (0.5, 1, 1), (1, 1, 1)),
                 'green': ((0, 0, 0), (0.5, 1, 1), (1, 0, 0)),
                 'blue':  ((0, 1, 1), (0.5, 1, 1), (1, 0, 0))}
        cmap = matplotlib.colors.LinearSegmentedColormap('cmap', cdict, 255)
        cs = axis.imshow(array, interpolation='nearest', extent=extent,
                                                     origin='lower', cmap=cmap)
        cs.set_clim(-max_val, max_val)
        axis.set_xlim(0.0, self.side_length)
        axis.set_ylim(0.0, self.side_length)
        toolbox_plotting.make_colorbar(axis, cs)
    def plot_population(self, axis, detail_level, sif_name):
        array = numpy.sign(self.get_rate_array(detail_level, sif_name))
        extent = [0.0, self.side_length]*2
        cdict = {'red':   ((0, 0, 0), (0.5, 1, 1), (1, 1, 1)),
                 'green': ((0, 0, 0), (0.5, 1, 1), (1, 0, 0)),
                 'blue':  ((0, 1, 1), (0.5, 1, 1), (1, 0, 0))}
        cmap = matplotlib.colors.LinearSegmentedColormap('cmap', cdict, 3)
        cs = axis.imshow(array, interpolation='nearest', extent=extent,
                                                     origin='lower', cmap=cmap)
        cs.set_clim(-1.0, 1.0)
        axis.set_xlim(0.0, self.side_length)
        axis.set_ylim(0.0, self.side_length)
    def calc_nearest_neighbor_distances(self, detail_level):
        rf = self.get_detail_level_results(detail_level)
        if rf.is_event_list_column_name('eNN_dist') and \
            rf.is_event_list_column_name('eNN_dist') and \
                rf.is_event_list_column_name('eNN_dist'):
            return rf
        cons_points = numpy.array([e.position() for e in rf.consumers()])
        prod_points = numpy.array([e.position() for e in rf.producers()])
        if self.wrapping:
            cons_points = wrap_points(cons_points, self.side_length, is_3d=self.is_3d)
            prod_points = wrap_points(prod_points, self.side_length, is_3d=self.is_3d)
        cons_tree = KDTree(cons_points)
        prod_tree = KDTree(prod_points)
        for e in rf.events:
            c_dist = cons_tree.query(e.position(), k=2)[0]
            p_dist= prod_tree.query(e.position(), k=2)[0]
            if (e.vars['mark'] == 2):
                e.vars['sNN_dist'] = c_dist[1]
                e.vars['oNN_dist'] = p_dist[0]
            else:
                e.vars['sNN_dist'] = p_dist[1]
                e.vars['oNN_dist'] = c_dist[0]
            e.vars['eNN_dist'] = min(e.vars['sNN_dist'], e.vars['oNN_dist'])
        rf.add_event_list_column_name('eNN_dist')
        rf.add_event_list_column_name('oNN_dist')
        rf.add_event_list_column_name('sNN_dist')
        #rf.eventList.update_text()
        rf.write()
        return rf
    def get_mean_NN_dist(self, detail_level, mean='amean', dist='oNN_dist'):
        rf = self.calc_nearest_neighbor_distances(detail_level)
        dists = [e.vars[dist] for e in rf.events]
        if not mean in ['amean', 'gmean', 'hmean']:
            toolbox_basic.error_message('toolbix_elmer.Simulation.get_mean_NN_dist()',
                        'mean not recognised: %s'%(mean))
        exec 'def mean(x): return %s(x)'%(mean)
        return mean(dists)
    def scatter_oNN_dist_vs_rate(self, axis, detail_level, sif_name, markersize=5):
        rf = self.get_combined_results(detail_level, sif_name)
        cons_rates = [-e.vars['rate'] for e in rf.consumers()]
        cons_dists = [e.vars['oNN_dist'] for e in rf.consumers()]
        axis.plot(cons_dists, cons_rates, '.', color='blue', markersize=markersize)
        prod_rates = [e.vars['rate'] for e in rf.producers()]
        prod_dists = [e.vars['oNN_dist'] for e in rf.producers()]
        axis.plot(prod_dists, prod_rates, '.', color='red', markersize=markersize)
        axis.set_xlabel(r'Interspecies N-N distance ($\mu$m)')
        axis.set_ylabel('Abs. metabolic rate '+r'(zmol $cell^{-1} ms^{-1}$)')
        axis.set_xlim(1, axis.get_xlim()[1])
    def scatter_oNN_dist_vs_concn(self, axis, detail_level, sif_name, markersize=5):
        rf = self.get_combined_results(detail_level, sif_name)
        cons_concns = [e.vars['amean_surf_concn'] for e in rf.consumers()]
        cons_dists = [e.vars['oNN_dist'] for e in rf.consumers()]
        axis.plot(cons_dists, cons_concns, '.', color='blue', markersize=markersize)
        prod_concns = [e.vars['amean_surf_concn'] for e in rf.producers()]
        prod_dists = [e.vars['oNN_dist'] for e in rf.producers()]
        axis.plot(prod_dists, prod_concns, '.', color='red', markersize=markersize)
        axis.set_xlabel('Interspecies nearest neighbour distance')
        axis.set_ylabel('Surface concentration')
        axis.set_xlim(1, axis.get_xlim()[1])
    def plot_kinetics(self, axis, detail_level, sif_name, maxs):
        consume_rate, produce_rate = \
                self.get_consume_produce_functions(detail_level, sif_name)
        p = list(numpy.linspace(0, maxs, num=1000))
        prod = [produce_rate(pval) for pval in p]
        cons = [-consume_rate(pval) for pval in p]
        axis.plot(p, prod, 'r-')
        axis.plot(p, cons, 'b-')
        #axis.set_xlabel(r'Product concentration ($\mu$M)')
        axis.set_xlabel(r'Hydrogen concentration ($\mu$M)')
        axis.set_ylabel('Metabolic rate '+r'(zmol $cell^{-1} ms^{-1}$)')
    def make_run_plot(self, detail_level, sif_name, maxP=None):
        fig = toolbox_plotting.ThesisFigure(double_column=True)
        axis = fig.add_subplot('A', 221)
        #self.plot_rate_array(axis, detail_level, sif_name)
        self.plot_population(axis, detail_level, sif_name)
        toolbox_plotting.empty_padding_axis(axis, "bottom")
        axis = fig.add_subplot('B', 222)
        toolbox_plotting.empty_padding_axis(axis, "bottom")
        if maxP == None:
            maxP = numpy.max(self.get_concn_array(detail_level, sif_name))
            maxP = 10**math.ceil(math.log10(maxP))
        self.plot_kinetics(axis, detail_level, sif_name, maxP)
        axis.text(8, 0.2, r'$q_{A}([H])$', color='r', va='center', ha='center')
        axis.text(8, 0.8, r'$-q_{B}([H])$', color='b', va='center', ha='center')
        analytic = AnalyticApproach()
        analytic.set_parameters(A=1, qmaxA=1, pmax=10, kA=10, qmaxB=5, pmin=0.04, kB=30)
        p_equal = analytic.calc_equal_concn()
        r_equal = analytic.production(p_equal)
        axis.plot([p_equal]*2, [0,r_equal+0.05],
                                    color='0.5', linestyle='-', zorder=-10)
        axis.text(p_equal, r_equal+0.05, '%.2f'%(p_equal),
                            color='0.5', va='bottom', ha='center', fontsize=8)
        axis.plot([0, p_equal+0.5], [r_equal]*2,
                                    color='0.5', linestyle='-', zorder=-10)
        axis.text(p_equal+0.6, r_equal, '%.2f'%(r_equal),
                             color='0.5', va='center', ha='left', fontsize=8)
        axis = fig.add_subplot('C', 223)
        cs = self.plot_concn_array(axis, detail_level, sif_name, plot_cs=False)
        cbar = toolbox_plotting.make_colorbar(axis, cs, side="bottom")
        #label = r'Product concentration ($\mu$M)'
        label = r'Hydrogen concentration ($\mu$M)'
        cbar.set_ticks([2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7])
        cbar.set_label(label)
        axis.set_xticklabels(['']*10)
        axis = fig.add_subplot('D', 224)
        toolbox_plotting.empty_padding_axis(axis, "bottom")
        self.calc_nearest_neighbor_distances(detail_level)
        self.scatter_oNN_dist_vs_rate(axis, detail_level, sif_name)
        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.96,
                                                    wspace=0.3, hspace=0.25)
        fig.process_subplots(label_pos=(0.0, 1.1))
        axis = fig.find_axis_from_label('C')
        axis.tick_params(bottom="off")
        fig.save(os.path.join(self.get_run_dir(detail_level, sif_name), 'run_plot.pdf'))
    def check_runs(self, rel_tol=1E-3):
        for dl_dir in toolbox_basic.subdir_list(self.root_dir, 'detail_level_*'):
            run_dirs = toolbox_basic.subdir_list(dl_dir)
            for run_dir in run_dirs:
                if os.path.basename(run_dir) == 'mesh':
                    continue
                cell_file = os.path.join(run_dir, 'cell_rates.xml')
                if os.path.isfile(cell_file):
                    cell_file = SimulationResultsFile(path=cell_file, read_only=True)
                    rel_diff = cell_file.get_relative_difference()
                    if rel_diff > rel_tol:
                        print('%s has rel_diff = %f'%(run_dir, rel_diff))
                else:
                    print('%s has no cell_rates.xml file'%(run_dir))
    def get_amean_concn(self, detail_level, sif_name):
        rf = self.get_run_results(detail_level, sif_name, read_only=True)
        return rf.get_amean_concn()
    def get_mean_surf_concn(self, detail_level, sif_name,
                                                cell_type='all', mean='amean'):
        rf = self.get_run_results(detail_level, sif_name, read_only=True)
        out = rf.get_mean_surf_concn(cell_type=cell_type, mean=mean)
        return out
    def get_sif_names(self, detail_level):
        dl_dir = self.get_detail_level_dir(detail_level)
        out = []
        for subdir in toolbox_basic.subdir_list(dl_dir):
            base_name = os.path.basename(subdir)
            sif_path = os.path.join(subdir, base_name+'.sif')
            if os.path.exists(sif_path):
                out.append(base_name)
        return out


class SimulationResultsFile(toolbox_results.ResultXMLfile):
    def __init__(self, path=None, read_only=False, header='mark,x,y'):
        toolbox_results.ResultXMLfile.__init__(self, path=path,
                                        root_name='elmer', read_only=read_only)
        self.simulation_root = self.find('./simulation')
        if self.simulation_root == None:
            self.simulation_root = xmlTree.SubElement(self.root, 'simulation')
        self.eventList = self.get_subresult('./eventList')
        if self.eventList == None:
            self.eventList = toolbox_result.ResultXMLfile(root_name='eventList')
            self.append_subresult(self.eventList)
        self.events = self.eventList.read_in_text()
        for line in self.events:
            line.__class__ = EventResult
    def get_event_list_column_names(self):
        return self.eventList.header.split(',')
    def set_space_parameters(self, wrapping=None, is_3d=None, side_length=None):
        space = self.find('./simulation/space')
        if space == None:
            sim = self.find('./simulation')
            space = xmlTree.SubElement(sim, 'space')
        if not wrapping == None: space.set('wrapping', str(wrapping))
        if not is_3d == None: space.set('is_3d', str(is_3d))
        if not side_length == None: space.set('side_length', str(side_length))
    def copy_event_list(self, simulation_results_file):
        self.remove_subresult(self.eventList)
        self.append_subresult(simulation_results_file.eventList)
        self.eventList = simulation_results_file.eventList
        self.events = simulation_results_file.events
    def set_event_list_header(self, header):
        self.eventList.header = header
        self.eventList.root.set('header', header)
    def add_event_list_column_name(self, column_name):
        if self.is_event_list_column_name(column_name): return
        self.eventList.add_column_name(column_name)
    def is_event_list_column_name(self, column_name):
        return (column_name in self.get_event_list_column_names())
    def set_up_population(self, detail_level, grid_length, wrapping, bio_array):
        for event in self.events:
            event.apply_detail_level(detail_level)
            event.setup_ranges(grid_length, wrapping)
            event.stamp_bio_array(bio_array)
    def apply_detail_level(self, detail_level):
        for event in self.events:
            event.apply_detail_level(detail_level)
    def setup_ranges(self, grid_length, wrapping):
        for event in self.events:
            event.setup_ranges(grid_length, wrapping)
    def stamp_bio_array(self, bio_array):
        for event in self.events:
            event.stamp_bio_array(bio_array)
    def consumers(self):
        return [e for e in self.events if e.vars['mark'] == 2]
    def producers(self):
        return [e for e in self.events if e.vars['mark'] == 3]
    def calc_rates_from_concn_array(self, concn_array, consume_rate,
                                                produce_rate, rate_array=None):
        for event in self.events:
            event.calc_rate_from_concn_array(concn_array, consume_rate,
                                          produce_rate, rate_array=rate_array)
            event.calc_amean_surf_concn(concn_array)
    def set_concn_rate_results(self, concn_array, rate_array):
        # This doesn't take account of detail_level!
        # See update_concn_rate_results()
        calculated_flux = \
                    numpy.sum(numpy.absolute(rate_array))/2
        rel_diff = abs(numpy.sum(rate_array))/calculated_flux
        max_concn = numpy.max(concn_array)
        amean_concn = amean(concn_array)
        min_concn = numpy.min(concn_array)
        concn_rate = self.find('./simulation/concn_rate')
        if concn_rate == None:
            sim = self.find('./simulation')
            concn_rate = xmlTree.SubElement(sim, 'concn_rate')
        concn_rate.set('calculated_flux', str(calculated_flux))
        concn_rate.set('rel_diff',  str(rel_diff))
        concn_rate.set('max_concn', str(max_concn))
        concn_rate.set('amean_concn', str(amean_concn))
        concn_rate.set('min_concn', str(min_concn))
    def update_concn_rate_results(self, detail_level):
        # This is a bit of a fudge: set_concn_rate_results() implicitly assumes
        # detail_level = 1
        production = numpy.sum([e.vars['rate'] for e in self.producers()])
        consumption = numpy.sum([e.vars['rate'] for e in self.consumers()])
        calculated_flux = (production - consumption)/2
        rel_diff = abs(production + consumption)/calculated_flux
        concn_rate = self.find('./simulation/concn_rate')
        concn_rate.set('calculated_flux', str(calculated_flux))
        concn_rate.set('rel_diff',  str(rel_diff))
    def get_amean_concn(self):
        return float(self.find('.simulation/concn_rate').attrib['amean_concn'])
    def get_calculated_flux(self):
        return float(self.find('.simulation/concn_rate').attrib['calculated_flux'])
    def get_relative_difference(self):
        return float(self.find('.simulation/concn_rate').attrib['rel_diff'])
    def get_mean_surf_concn(self, cell_type='all', mean='amean'):
        if cell_type == 'all':
            events = self.events
        elif cell_type == 'consumers':
            events = self.consumers()
        elif cell_type == 'producers':
            events = self.producers()
        if not mean in ['amean', 'gmean', 'hmean']:
            toolbox_basic.error_message('toolbix_elmer.Simulation.get_mean_NN_dist()',
                        'mean not recognised: %s'%(mean))
        exec 'def mean(x): return %s(x)'%(mean)
        return mean([e.vars['amean_surf_concn'] for e in events])





class EventResult(toolbox_results.SingleCSVline):
    def __init__(self, header, text):
        toolbox_results.SingleCSVline.__init__(self, header, text)
    def apply_detail_level(self, detail_level):
        # If detail_level is odd:
        if (detail_level%2 == 1):
            diff = int((detail_level-1)/2)
            i_cen = int(self.vars['x'] * detail_level)
            self.vars['x'] = (i_cen+0.5)/detail_level
            self.vars['i_min'] = i_cen - diff
            self.vars['i_max'] = i_cen + diff + 1
            j_cen = int(self.vars['y'] * detail_level)
            self.vars['y'] = (j_cen+0.5)/detail_level
            self.vars['j_min'] = j_cen - diff
            self.vars['j_max'] = j_cen + diff + 1
            if 'z' in self.vars.keys():
                k_cen = int(self.vars['z'] * detail_level)
                self.vars['z'] = (k_cen+0.5)/detail_level
                self.vars['k_min'] = k_cen - diff
                self.vars['k_max'] = k_cen + diff + 1
        # If detail_level is even:
        else:
            diff = int(detail_level/2)
            i_cen = int(round(self.vars['x'] * detail_level))
            self.vars['x'] = i_cen/detail_level
            self.vars['i_min'] = i_cen - diff
            self.vars['i_max'] = i_cen + diff
            j_cen = int(round(self.vars['y'] * detail_level))
            self.vars['y'] = j_cen/detail_level
            self.vars['j_min'] = j_cen - diff
            self.vars['j_max'] = j_cen + diff
            if 'z' in self.vars.keys():
                k_cen = int(round(self.vars['z'] * detail_level))
                self.vars['z'] = k_cen/detail_level
                self.vars['k_min'] = k_cen - diff
                self.vars['k_max'] = k_cen + diff
    def setup_ranges(self, grid_length, wrapping):
        # Take care of any edge effects:
        i_range = range(self.vars['i_min'], self.vars['i_max'])
        j_range = range(self.vars['j_min'], self.vars['j_max'])
        if 'z' in self.vars.keys():
            k_range = range(self.vars['k_min'], self.vars['k_max'])
        if wrapping:
            self.i_range = [i%grid_length for i in i_range]
            self.j_range = [j%grid_length for j in j_range]
            if 'z' in self.vars.keys():
                self.k_range = [k%grid_length for k in k_range]
        else:
            self.i_range = [i for i in i_range if i >= 0 and i <= grid_length]
            self.j_range = [j for j in j_range if j >= 0 and j <= grid_length]
            if 'z' in self.vars.keys():
                self.k_range = [k for k in k_range if k>=0 and k<=grid_length]
    def stamp_bio_array(self, bio_array):
        for i in self.i_range:
            for j in self.j_range:
                if 'z' in self.vars.keys():
                    for k in self.k_range:
                        bio_array[i][j][k] = self.vars['mark']
                else:
                    bio_array[i][j] = self.vars['mark']
    def calc_rate_from_concn_array(self, concn_array, consume_rate,
                                        produce_rate, rate_array=None):
        self.vars['rate'] = 0.0
        if self.vars['mark'] == 2: kinetic_rate = consume_rate
        else:                      kinetic_rate = produce_rate
        counter = 0
        for (i, j) in [(i, j) for i in self.i_range for j in self.j_range]:
            concns = [concn_array[I][J] for I in [i, i+1] for J in [j, j+1]]
            rates = [kinetic_rate(concn) for concn in concns]
            mean_rate = numpy.mean(rates)
            if not rate_array == None:
                rate_array[i][j] = mean_rate
            self.vars['rate'] += mean_rate
            counter += 1
        self.vars['rate'] /= counter
        return self.vars['rate']
    def calc_amean_surf_concn(self, concn_array):
        surface_nodes = [(i, self.j_range[0]) for i in self.i_range] + \
                     [(i, self.j_range[-1]+1) for i in self.i_range] + \
                     [(self.i_range[0], j) for j in self.j_range[1:]] + \
           [(self.i_range[-1]+1, j) for j in self.j_range+[self.j_range[-1]+1]]
        concns = [concn_array[i][j] for (i, j) in surface_nodes]
        self.vars['amean_surf_concn'] = numpy.mean(concns)
        return self.vars['amean_surf_concn']
    def position(self):
        if 'z' in self.vars.keys():
            return (self.vars['x'], self.vars['y'], self.vars['z'])
        else:
            return (self.vars['x'], self.vars['y'])


class AnalyticApproach:
    def __init__(self):
        self.A = 1
        self.qmaxA = 1.0
        self.pmax = 10.0
        self.kA = 1.0
        self.B = 1
        self.qmaxB = 1.0
        self.pmin = 0.1
        self.kB = 1.0
    def set_parameters(self, A=None, qmaxA=None, pmax=None, kA=None,
                             B=None, qmaxB=None, pmin=None, kB=None):
        self.A = self.A if A == None else A
        self.qmaxA = self.qmaxA if qmaxA == None else qmaxA
        self.pmax = self.pmax if pmax == None else pmax
        self.kA = self.kA if kA == None else kA
        self.B = self.B if B == None else B
        self.qmaxB = self.qmaxB if qmaxB == None else qmaxB
        self.pmin = self.pmin if pmin == None else pmin
        self.kB = self.kB if kB == None else kB
    def production(self, p):
        return self.A*self.qmaxA*(self.pmax-p)/(self.pmax+self.kA+p)
    def calc_equal_concn(self):
        qmaxAA, qmaxBB = self.qmaxA*self.A, self.qmaxB*self.B
        q2 = qmaxAA + qmaxBB
        q1 = qmaxAA*(self.kB + self.pmin - self.pmax) \
                + qmaxBB*(self.kA + self.pmax - self.pmin)
        q0 = - qmaxBB*self.kA*self.pmin \
                - qmaxAA*self.kB*self.pmax - q2*self.pmax*self.pmin
        roots = numpy.roots([q2, q1, q0])
        p = max(roots)
        return p
    def calc_equal_concn_rate(self):
        p = self.calc_equal_concn()
        return self.production(p)
    def sensitivity_analysis(self, cv=0.1, return_rate=False, return_diffs=True):
        params = (self.A, self.qmaxA, self.pmax, self.kA,
                  self.B, self.qmaxB, self.pmin, self.kB)
        if return_rate:
            norm_val = self.calc_equal_concn_rate()
        else:
            norm_val = self.calc_equal_concn()
        max_val, min_val = norm_val, norm_val
        cv_range = [(1-cv), 1, (1+cv)]
        for a in cv_range:
            for qa in cv_range:
                for px in cv_range:
                    for ka in cv_range:
                        for b in cv_range:
                            for qb in cv_range:
                                for pn in cv_range:
                                    for kb in cv_range:
                                        self.set_parameters(A=a*params[0],
                                                            qmaxA=qa*params[1],
                                                            pmax=px*params[2],
                                                            kA=ka*params[3],
                                                            B=b*params[4],
                                                            qmaxB=qb*params[5],
                                                            pmin=pn*params[6],
                                                            kB=kb*params[7])
                                        if return_rate:
                                            val = self.calc_equal_concn_rate()
                                        else:
                                            val = self.calc_equal_concn()
                                        max_val = max(max_val, val)
                                        min_val = min(min_val, val)
        self.set_parameters(A=params[0], qmaxA=params[1], pmax=params[2],
                            kA=params[3], B=params[4], qmaxB=params[5],
                            pmin=params[6], kB=params[7])
        if return_diffs:
            minus_diff = norm_val - min_val
            plus_diff = max_val - norm_val
            return minus_diff, plus_diff
        else:
            return min_val, max_val


class SimCollection:
    def __init__(self, simulation):
        if isinstance(simulation, list):
            self.simulations = simulation
            simulation = self.simulations[0]
        else:
            self.simulations = [simulation]
        self.wrapping = simulation.wrapping
        self.side_length = simulation.side_length
        self.is_3d = simulation.is_3d
        self.consumers = simulation.consumers
        self.producers = simulation.producers
        self.random_seed = simulation.random_seed
    def add_if_belongs(self, simulation, diffs_allowed=['random_seed']):
        comparitor = self.simulations[0]
        belongs = True
        if not simulation.wrapping == self.wrapping and \
                                    not 'wrapping' in diffs_allowed:
            belongs = False
        if not simulation.side_length == self.side_length and \
                                    not 'side_length' in diffs_allowed:
            belongs = False
        if not simulation.is_3d == self.is_3d and \
                                    not 'is_3d' in diffs_allowed:
            belongs = False
        if not simulation.consumers == self.consumers and \
                                    not 'consumers' in diffs_allowed:
            belongs = False
        if not simulation.producers == self.producers and \
                                    not 'producers' in diffs_allowed:
            belongs = False
        if not simulation.random_seed == self.random_seed and \
                                    not 'random_seed' in diffs_allowed:
            belongs = False
        if belongs:
            self.simulations.append(simulation)
        return belongs
    def get_calculated_fluxes(self, detail_level, sif_name):
        out = []
        for sim in self.simulations:
            rf = sim.get_run_results(detail_level, sif_name)
            if detail_level > 1:
                rf.update_concn_rate_results(detail_level)
            out.append(rf.get_calculated_flux())
        return out
    def get_amean_concns(self, detail_level, sif_name):
        return \
            [sim.get_run_results(detail_level, sif_name).get_amean_concn() \
                                                for sim in self.simulations]
    def estimates_from_concn(self, detail_level, sif_name, D, pmin,
                                                    dist_mean='amean'):
        sides = 6 if self.is_3d else 4
        estimates = []
        for sim in self.simulations:
            p = sim.get_amean_concn(detail_level, sif_name)
            d = sim.get_mean_NN_dist(detail_level, mean=dist_mean)
            estimates.append(D*self.producers*sides*(p-pmin)/d)
        return estimates
    def estimates_from_surf_concn(self, detail_level, sif_name, D,
                                                    dist_mean='amean'):
        sides = 6 if self.is_3d else 4
        estimates = []
        for sim in self.simulations:
            pmin = sim.get_mean_surf_concn(detail_level, sif_name,
                                        cell_type='consumers', mean='amean')
            pmax = sim.get_mean_surf_concn(detail_level, sif_name,
                                        cell_type='producers', mean='amean')
            d = sim.get_mean_NN_dist(detail_level, mean=dist_mean)
            estimates.append(D*self.producers*sides*(pmax-pmin)/d)
        return estimates


def find_sim_collections(results_dir, diffs_allowed=['random_seed']):
    sim_collections = []
    for sim in get_replicate_simulations(results_dir):
        sim_collection = None
        for sc in sim_collections:
            if sc.add_if_belongs(sim):
                sim_collection = sc
                break
        if sim_collection == None:
            sim_collection = SimCollection(sim)
            sim_collections.append(sim_collection)
    return sim_collections


# points should be given as a 2xn numpy.array
# side_length should be a positive real number (usually an integer)
def wrap_points(points, side_length, num_wraps=1, is_3d=False):
    diffs = [i*side_length for i in range(num_wraps+1)]
    diffs.extend([-d for d in diffs[1:]])
    new_points = []
    for point in points:
        for x_diff in diffs:
            for y_diff in diffs:
                if is_3d:
                    for z_diff in diffs:
                        new_points.append(point + (x_diff, y_diff, z_diff))
                else:
                    new_points.append(point + (x_diff, y_diff))
    return numpy.array(new_points)


def every_cell_spatial_stats(root_dir, detail_level=1):
    results_file = get_detail_level_results(root_dir, detail_level=detail_level)
    results_file.get_space_parameters()
    points = numpy.array([e.position() for e in results_file.events])
    if results_file.wrapping:
        points = wrap_points(points, results_file.side_length,
                                                    is_3d=results_file.is_3d)
    triangulation = Delaunay(points)
    indices, indptr = triangulation.vertex_neighbor_vertices
    for event in results_file.events:
        e_point = numpy.array(event.position())
        row_number = toolbox_basic.find_index_of_row_in_array(points, e_point)
        event.vars['eDT_nnbs'] = indices[row_number+1] - indices[row_number] +1
        neighbor_indices = indptr[indices[row_number]:indices[row_number+1]]
        nb_points = [points[i] for i in neighbor_indices]


def get_replicate_results(replicates_dir):
    replicates_dir = toolbox_basic.check_path(replicates_dir)
    results_path = os.path.join(replicates_dir, 'results.xml')
    results_file = ReplicatesFile(path=results_path)
    return results_file


def get_replicate_simulations(replicates_dir):
    return [Simulation(d) for d in toolbox_basic.subdir_list(replicates_dir)]


def setup_replicate_results(replicates_dir):
    results_file = get_replicate_results(replicates_dir)
    replicate_simulations = get_replicate_simulations(replicates_dir)
    wrapping = replicate_simulations[0].wrapping
    is_3d=replicate_simulations[0].is_3d
    side_length=replicate_simulations[0].side_length
    results_file.set_space_parameters(wrapping=wrapping, is_3d=is_3d, side_length=side_length)
    for sim in replicate_simulations:
        if not ((sim.wrapping == wrapping) and
                (sim.is_3d == is_3d) and
                (sim.side_length == side_length)):
            toolbox_basic.error_message('toolbox_elmer.get_replicates_results():'+
                'Replicates have different space parameters', replicates_dir)
            exit()
    results_file.write()
    return results_file


def get_replicate_results_basics(replicates_dir,
                                    detail_level=1, sif_name='elmer_standard'):
    replicate_simulations = get_replicate_simulations(replicates_dir)
    results_file = setup_replicate_results(replicates_dir)
    dl_resuls = results_file.get_detail_level_results(detail_level=detail_level)
    sif_results = results_file.get_solver_input_file_results(
                                  detail_level=detail_level, sif_name=sif_name)
    #for sim in replicate_simulations:

    results_file.write()



def hydrogen_logo(axis, bottom_left=(0.9, 0.9), height_width=0.1):
    color = '0.5'
    circle = toolbox_schematic.Circle()
    circle.set_defaults(edgecolor='none', facecolor=color, transform=True)
    radius = 0.2*height_width
    center_A = (bottom_left[0] + radius, bottom_left[1] + height_width - radius)
    center_B = (bottom_left[0] + height_width - radius, bottom_left[1] + radius)
    circle.set_points(center_A, radius)
    circle.draw(axis)
    circle.set_points(center_B, radius)
    circle.draw(axis)
    print center_A, center_B
    axis.plot([center_A[0], center_B[0]], [center_A[1], center_B[1]], color, linestyle='-', transform=axis.transAxes)
    #center = ((bottom_left[0]+top_right[0])/2, (bottom_left[1]+top_right[1])/2)




'''
def get_combined_results(root_dir, detail_level=1, sif_name='elmer_standard'):
    dl_results = get_detail_level_results(root_dir, detail_level=detail_level)
    run_results = get_run_results(root_dir,
                                  detail_level=detail_level, sif_name=sif_name)
    for e in run_results.events:
        for o in dl_results.events:
            if o.position() == e.position():
                for attrib in o.vars.keys():
                    e.vars[attrib] = o.vars[attrib]
    run_results.set_read_only()
    return run_results
'''

'''
def get_rate_array(root_dir, detail_level=1, sif_name='elmer_standard'):
    run_dir = get_run_dir(root_dir,detail_level=detail_level,sif_name=sif_name)
    array_path = os.path.join(run_dir, 'rate_array.npy')
    array_path = toolbox_basic.check_path(array_path)
    array = numpy.load(array_path)
    return array
'''

'''
def calc_nearest_neighbor_distances(root_dir, detail_level=1):
    simulation = Simulation(root_dir)
    results_file = simulation.get_detail_level_results(detail_level)
    #if ('eNN_dist' in results_file.get_eventList_column_names()): return
    if results_file.is_event_list_column_name('eNN_dist'): return
    #results_file.set_space_parameters()
    cons_points = numpy.array([e.position() for e in results_file.consumers()])
    prod_points = numpy.array([e.position() for e in results_file.producers()])
    if simulation.wrapping:
        cons_points = wrap_points(cons_points, simulation.side_length,
                                                    is_3d=simulation.is_3d)
        prod_points = wrap_points(prod_points, simulation.side_length,
                                                    is_3d=simulation.is_3d)
    cons_tree = KDTree(cons_points)
    prod_tree = KDTree(prod_points)
    for e in results_file.events:
        c_dist, id = cons_tree.query(e.position())
        p_dist, id = prod_tree.query(e.position())
        e.vars['sNN_dist'] = c_dist if (e.vars['mark'] == 2) else p_dist
        e.vars['oNN_dist'] = c_dist if (e.vars['mark'] == 3) else p_dist
        e.vars['eNN_dist'] = min(c_dist, p_dist)
    results_file.add_event_list_column_name('eNN_dist')
    results_file.add_event_list_column_name('oNN_dist')
    results_file.add_event_list_column_name('sNN_dist')
    results_file.eventList.update_text()
    results_file.write()
'''
