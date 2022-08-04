#!/usr/bin/python
from __future__ import division
import math
import numpy
import os
import toolbox_results
import toolbox_schematic

#This has currently been modified to do results for 10 days instead of the usual 100
# By convention, in the main set of results the first character of the 
# protocol/directory name refers to the segregation strategy, the second to the 
# investment in repair, and the third & fourth to the damage accumulation rate
# a: asymmetry, m: mid-symmetry, s: symmetry
# n: no repair, r: fixed repair, o: optimal repair
# accumulation rates are multiplied by 100, so 15 refers to 0.15 h-1
def process_default_results(results_path, attribute):
    attributes = {'header':'directory,mean,std',
                  'starting_time':'240', 'name':attribute}
    results_output = toolbox_results.ResultsOutput(path=results_path)
    results_set = toolbox_results.ResultSet(results_output, attributes)
    for result in results_set.members:
        dir = result.vars['directory']
        result.vars['segregation'] = dir[0]
        result.vars['repair'] = dir[1]
        result.vars['accumulation'] = 0.01*float(dir[2:4])
        result.vars['mean'] = float(result.vars['mean'])
        result.vars['std'] = float(result.vars['std'])
        # Convert from g L-1 to mg L-1
        if attribute == 'glucose' or attribute == 'totalBiomass':
            result.vars['mean'] *= 1000
            result.vars['std'] *= 1000
    return results_set


def process_erjavec_results(results_path):
    attributes = {'header':'directory,mean,std',
                  'starting_time':'240', 'name':'specific growth rate'}
    results_output = toolbox_results.ResultsOutput(path=results_path)
    results_set = toolbox_results.ResultSet(results_output, attributes)
    for result in results_set.members:
        dir = result.vars['directory']
        result.vars['segregation'] = dir[0]
        result.vars['repair'] = 'n'
        result.vars['accumulation'] = 0.1*float(dir[1:3])
        result.vars['mean'] = float(result.vars['mean'])
        result.vars['std'] = float(result.vars['std'])
    return results_set


# This is for the set of results where investment in repair, beta, was varied 
def process_beta_results(results_path):
    attributes = {'header':'directory,mean,std',
                  'starting_time':'240', 'name':'specific growth rate'}
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    segregation = []
    for result in result_set.members:
        dir = result.vars['directory']
        # 'a', 'm', or 's'
        result.vars['segregation'] = dir[0]
        if not dir[0] in segregation: segregation.append(dir[0])
        # 'repair' is 'n' to ensure a solid line in the plot
        result.vars['repair'] = 'n'
        # 'accumulation' is actually the value of beta
        result.vars['accumulation'] = 0.01 * float(dir[2:4])
        # gGrowth rate mean & standard deviation
        result.vars['mean'] = float(result.vars['mean'])
        result.vars['std'] = float(result.vars['std'])
        result.vars['optimal'] = False
    for seg in segregation:
        seg_set = [r for r in result_set.members if r.vars['segregation']==seg]
        means = [r.vars['mean'] for r in seg_set]
        seg_set[means.index(max(means))].vars['optimal'] = True
    return result_set


# Take the glucose and totalBiomass concentrations and calculate the yield
# Yield is calculated using the formula 
#        Y = biomass / (inflow_concentration - glucose_concentration)
# Standard deviation is calculated using the formula 
#        CV(x/y) = [CV(x)^2 + CV(y)^2]^0.5
def process_yield_results(results_path, inflow=3.24):
    glucose_results = process_default_results(results_path, 'glucose')
    biomass_results = process_default_results(results_path, 'totalBiomass')
    attributes = {'header':'directory,mean,std',
                  'starting_time':'240', 'name':'yield'}
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for g_result in glucose_results.members:
        dirname = g_result.vars['directory']
        mean_g = g_result.vars['mean']
        std_g = g_result.vars['std']
        b_result = biomass_results.find_single_result('directory', dirname)
        mean_b = b_result.vars['mean']
        std_b = b_result.vars['std']
        y_result = toolbox_results.SingleResult()
        y_result.vars['directory'] = dirname
        y_result.vars['accumulation'] = g_result.vars['accumulation']
        y_result.vars['segregation'] = g_result.vars['segregation']
        y_result.vars['repair'] = g_result.vars['repair']
        if mean_b == 0.0:
            y_result.vars['mean'] = 'nan'
            y_result.vars['std'] = 'nan'
        else:
            mean_y = mean_b / (inflow - mean_g)
            cv_g = (std_g / mean_g)
            cv_b = (std_b / mean_b)
            y_result.vars['mean'] = mean_y
            y_result.vars['std'] = mean_y * math.sqrt(cv_g**2 + cv_b**2) 
        result_set.members.append(y_result)
    return result_set    


def process_optima_results(results_path, attribute):
    attributes = {'header':'directory,mean,std',
                  'starting_time':'240', 'name':'optimum '+attribute}
    results_output = toolbox_results.ResultsOutput(path=results_path)
    result_set = toolbox_results.ResultSet(results_output, attributes)
    for result in result_set.members:
        dir = result.vars['directory']
        result.vars['segregation'] = dir[0]
        result.vars['repair'] = 'n'
        # This is actually the value of beta
        result.vars['accumulation'] = 0.01 * float(dir[5:7])
        result.vars['mean'] = 0.01 * float(dir[2:4])
        result.vars['std'] = 0.005
    return result_set


def plot_results(axis, results_set, all_acc=True,
                         colors={'a':'red','m':'green','s':'blue'},
                         linestyles={'n':'-','r':'--','o':'-.'}, error=True):
    for segregation, color in colors.iteritems():
        for repair, linestyle in linestyles.iteritems():
            #print segregation, repair
            data = [r for r in results_set.members if 
                                    r.vars['segregation'] == segregation and
                                    r.vars['repair'] == repair]
            data.sort(key=lambda r: r.vars['accumulation'])
            if all_acc:
                x = [r.vars['accumulation'] for r in data]
            else:
                x = [0.00, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
                data = [r for r in data if x.count(r.vars['accumulation'])>0 ]
            y = [r.vars['mean'] for r in data]
            yerr = [r.vars['std'] for r in data]
            if y.count('nan') > 0:
                stop_index = y.index('nan')
            elif x.count(0.26) > 0:
                stop_index = x.index(0.26)
            else:
                stop_index = None
            #print x[:stop_index], y[:stop_index]
            axis.plot(x[:stop_index], y[:stop_index],
                            color=color,linestyle=linestyle,zorder=5)
            if error:
                axis.errorbar(x[:stop_index], y[:stop_index], 
                            yerr=yerr[:stop_index], ecolor='grey', 
                            fmt=None, elinewidth=0.5, capsize=1, zorder=1)


def zoom_inset(axis, results_set, x_ticks, y_ticks, all_acc=True,
            source={'x_min':0.0, 'x_max':0.25, 'y_min':0.0, 'y_max':0.6},
            target={'x_left':0.0, 'x_right':0.25, 'y_top':0.6, 'y_bottom':0.0},
            colors={'a':'red','m':'green','s':'blue'},
            linestyles={'n':'-','r':'--','o':'-.'}, error=True):
    source_width = source['x_max'] - source['x_min']
    source_height = source['y_max'] - source['y_min']
    target_width = target['x_right'] - target['x_left']
    target_height = target['y_top'] - target['y_bottom']
    tick_font_size = 8
    for segregation, color in colors.iteritems():
        for repair, linestyle in linestyles.iteritems():
            interpolate_left, interpolate_right = False, False
            data = [r for r in results_set.members if 
                                    r.vars['segregation'] == segregation and
                                    r.vars['repair'] == repair]
            data.sort(key=lambda r: r.vars['accumulation'])
            if all_acc:
                x = [r.vars['accumulation'] for r in data]
            else:
                x = [0.00, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
                data = [r for r in data if r.vars['accumulation'] in x]
            # Check if we need to interpolate on the left
            if not source['x_min'] in x:
                interpolate_left = True
                result_left = [d for d in data if 
                                d.vars['accumulation'] < source['x_min']][-1]
                result_right = [d for d in data if 
                                d.vars['accumulation'] > source['x_min']][0]
                acc_left = result_left.vars['accumulation']
                mean_left = result_left.vars['mean']
                acc_right = result_right.vars['accumulation']
                mean_right = result_right.vars['mean']
                factor = (source['x_min'] - acc_left) / (acc_right - acc_left)
                interpolated_mean = mean_left + (mean_right - mean_left)*factor
                interpolate = toolbox_results.SingleResult()
                interpolate.vars['mean'] = interpolated_mean
                interpolate.vars['accumulation'] = source['x_min']
                data.insert(0,interpolate) 
            data = [r for r in data if 
                            r.vars['accumulation'] >= source['x_min']]
            data.sort(key=lambda r: r.vars['accumulation'])
            x = [r.vars['accumulation'] for r in data]
            # Check if we need to interpolate on the right
            if not source['x_max'] in x:
                interpolate_right = True
                result_left = [d for d in data if 
                                d.vars['accumulation'] < source['x_max']][-1]
                result_right = [d for d in data if 
                                d.vars['accumulation'] > source['x_max']][0]
                acc_left = result_left.vars['accumulation']
                mean_left = result_left.vars['mean']
                acc_right = result_right.vars['accumulation']
                mean_right = result_right.vars['mean']
                factor = (source['x_max'] - acc_left) / (acc_right - acc_left)
                interpolated_mean = mean_left + (mean_right - mean_left)*factor
                interpolate = toolbox_results.SingleResult()
                interpolate.vars['mean'] = interpolated_mean
                interpolate.vars['accumulation'] = source['x_max']
                data.append(interpolate)
            data = [r for r in data if 
                            r.vars['accumulation'] <= source['x_max']]
            data.sort(key=lambda r: r.vars['accumulation'])
            x = [r.vars['accumulation'] for r in data]
            # Scale the data
            for r in data:
                acc = r.vars['accumulation'] - source['x_min']
                acc *= target_width / source_width
                r.vars['accumulation'] = acc + target['x_left']
                mean = r.vars['mean'] - source['y_min']
                mean *= target_height / source_height
                r.vars['mean'] = mean + target['y_bottom']
                if 'std' in r.vars.keys():
                    std = r.vars['std'] * target_height / source_height
            x = [r.vars['accumulation'] for r in data]
            y = [r.vars['mean'] for r in data]
            if y.count('nan') > 0:
                stop_index = y.index('nan')
            elif x.count(0.26) > 0:
                stop_index = x.index(0.26)
            else:
                stop_index = None
            axis.plot(x[:stop_index], y[:stop_index],
                            color=color,linestyle=linestyle,zorder=5)
            start_index = 1 if interpolate_left else 0
            if interpolate_right:
                if stop_index == None: stop_index = -1
                else:                  stop_index -= 1
            yerr = [r.vars['std'] for r in data[start_index:stop_index]]
            if error:
                axis.errorbar(x[start_index:stop_index], y[start_index:stop_index], 
                            yerr=yerr, ecolor='grey', 
                            fmt=None, elinewidth=0.5, capsize=1, zorder=1)
            axis.plot([source['x_min'], source['x_max'], source['x_max'], source['x_min'], source['x_min']],
                        [source['y_min'], source['y_min'], source['y_max'], source['y_max'], source['y_min'], ],
                        '0.5', linewidth=1, zorder=-100)
            ax_left = target['x_left'] - 0.02*target_width
            ax_right = target['x_right'] + 0.02*target_width
            ax_bottom = target['y_bottom'] - 0.02*target_height
            ax_top = target['y_top'] + 0.02*target_height
            axis.plot([ax_left, ax_right, ax_right, ax_left, ax_left],
                        [ax_bottom, ax_bottom, ax_top, ax_top, ax_bottom],
                        '0.5', linewidth=1)
            tick_length = 0.005
            for tick in y_ticks:
                y = target['y_bottom'] + (tick - source['y_min'])*target_height/source_height
                axis.plot([ax_left, ax_left - tick_length], [y, y], '0.5', linewidth=1)
                axis.text(ax_left - 2*tick_length, y, str(tick), color='0.5',
                            ha='right', va='center', fontsize=tick_font_size)
            y2xscale = axis.get_ylim()[1] / axis.get_xlim()[1]
            tick_length *= y2xscale
            for tick in x_ticks:
                x = target['x_left'] + (tick - source['x_min'])*target_width / source_width
                axis.plot([x, x], [ax_bottom, ax_bottom - tick_length], '0.5', linewidth=1)
                axis.text(x, ax_bottom - 2*tick_length, str(tick), color='0.5',
                                ha='center', va='top', fontsize=tick_font_size)


def scatter_population(axis, species_output, x_attribute, y_attribute,
                    color='k', style='o', markersize=2):
    x_vals, y_vals = [], []
    biomass_names = ['activeBiomassGrowth', 'activeBiomassRepair', 'inactiveBiomassGrowth', 'inactiveBiomassRepair']
    for agent in species_output.members:
        if x_attribute == 'totalBiomass':
            x_vals.append(agent.get_total_biomass(biomass_names))
        elif x_attribute == 'specificGrowthRate':
            x_vals.append(agent.get_specific_growth_rate(biomass_names))
        else:
            x_vals.append(float(agent.vars[x_attribute]))
        if y_attribute == 'totalBiomass':
            y_vals.append(agent.get_total_biomass(biomass_names))
        elif y_attribute == 'specificGrowthRate':
            y_vals.append(agent.get_specific_growth_rate(biomass_names))
        else:
            y_vals.append(float(agent.vars[y_attribute]))
    axis.plot(x_vals, y_vals, style, color=color, markeredgecolor='none', markersize=markersize)


def draw_cell(axis, x1, x2, y1, r, arrow=0.003,
                        y2xscale=1, growth=False, repair=1, toxic=1):
    black = '0.5'
    cell = toolbox_schematic.Cell()
    cell.set_defaults(y2xscale=y2xscale)
    cell.set_points([x1, y1], [x2, y1], r)
    cell.draw(axis, edgecolor=black)
    # S -> Pact
    x_left, x_right, y = x1 - 1.5*r, x1 + r/2, y1 + y2xscale * r/5
    axis.plot([x_left, x_right], [y]*2, black)
    axis.plot([x_right - arrow, x_right, x_right - arrow],
            [y + y2xscale*arrow, y, y - y2xscale*arrow], black)
    # Pact -> Pdam
    x, y_top, y_bottom = x1 + r, y1, y1 - y2xscale * r/2
    axis.plot([x]*2, [y_top, y_bottom], black)
    axis.plot([x - arrow, x, x + arrow],
                [y_bottom + y2xscale*arrow, y_bottom, y_bottom + y2xscale*arrow], black)
    # Pdam -> Pact
    x = x + 2.5*arrow
    axis.plot([x]*2, [y_top, y_bottom], black)
    axis.plot([x - arrow, x, x + arrow],
                        [y_top - y2xscale*arrow, y_top, y_top - y2xscale*arrow], black)
    # Growth
    x_left, x_right = x1, x + 2* arrow
    y_top, y_bottom = y1 + 0.6 * y2xscale * r, y1 + 0.4 * y2xscale * r
    c = 'blue' if growth else black
    axis.plot([x_right, x_left, x_left], [y_top, y_top, y_bottom], color=c)
    # Repair
    x_left, x_right = x + arrow, x + 2* arrow
    y_top, y_bottom = y, y1 - 0.4 * y2xscale * r
    if repair == 2:   c = 'green'
    elif repair == 1: c = black
    else:             c = 'white'
    axis.plot([x_right, x_right, x_left], [y_top, y_bottom, y_bottom], color=c)
    # Toxicity
    x_right, x_left = x1 + 0.6 * r, x1
    y_bottom, y_top = y1 - 0.6 * y2xscale * r, y1
    if toxic:
        axis.plot([x_right, x_left, x_left, x_left - arrow/2, x_left + arrow/2],
            [y_bottom, y_bottom, y_top, y_top, y_top], 'r-')


def draw_chemostat(axis, x1, y1, x2, y2, color='0.5'):
    # Draw the inflow arrow
    x_temp = (7.5*x1 + 2.5*x2)/10
    y_temp = (3*y1 + 7*y2)/10
    axis.plot([x1, x_temp, x_temp], [y2, y2, y_temp], color)
    axis.plot([x_temp], [y_temp], 'v', markeredgecolor='none', markerfacecolor=color)
    # Draw the outflow arrow
    x_temp = (2.5*x1 + 7.5*x2)/10
    axis.plot([x_temp, x_temp, x2], [y_temp, y2, y2], color)
    axis.plot([x2], [y2], '>', markeredgecolor='none', markerfacecolor=color)
    # Draw the 3-sided box
    x1_temp = (9*x1 + 1*x2)/10
    x2_temp = (1*x1 + 9*x2)/10
    y_temp = (2*y1 + 8*y2)/10
    axis.plot([x1_temp, x1_temp, x2_temp, x2_temp], [y_temp, y1, y1, y_temp], color)
    # Draw the stirrer
    x_temp = (x1 + x2)/2
    y_temp = (6.5*y1 + 3.5*y2)/10
    x_vals, y_vals = [x_temp], [y2]
    x_const, y_const = (x2-x1)/5, (y2-y1)/5
    for angle in numpy.linspace(-0.5*math.pi, 1.5*math.pi, 100):
        sin, cos = math.sin(angle), math.cos(angle)
        x_vals.append(x_temp + (x_const*cos/(1+sin**2)))
        y_vals.append(y_temp + (y_const*sin*cos/(1+sin**2)))
    axis.plot(x_vals, y_vals, color)


def draw_const_env(axis, x1, y1, x2, y2, color='0.5', y2xscale=0.6/0.25):
    rscale, arrow = 0.75, 0.003
    radius, height = (x2 - x1)/(10), y2 - y1
    arrow = height/20
    # Draw the newly-divided cell
    x = x1 + 3*radius
    y_bottom, y_top = y1 + 0.6*height + y2xscale*radius, y1 + height - y2xscale*radius
    cell = toolbox_schematic.Cell()
    cell.set_defaults(edgecolor=color, y2xscale=y2xscale)
    cell.set_points([x, y_top], [x, y_bottom], radius*rscale)
    cell.draw(axis)
    # Draw the layer of cells
    y_bottom, y_top = y1 + 0.2*height + y2xscale*radius, y1 + 0.6*height - y2xscale*radius
    for i in range(5):
        x = x1 + (2 * i + 1) * radius
        cell.set_points([x, y_top], [x, y_bottom], radius*rscale)
        cell.draw(axis)
        axis.plot([x, x], [y1, y_bottom - y2xscale*radius - arrow], color=color)
        axis.plot([x], [y_bottom - y2xscale*radius - arrow], '^',
                    markerfacecolor=color, markeredgecolor='none', markersize=4)
    # The substratum
    axis.plot([x1, x2], [y1]*2, color=color)
    # Show the new cell replacing the displaced cell
    x_cen, y_top = (x1+x2)/2, y1 + 0.7*height
    x_vals, y_vals = [], []
    for angle in numpy.linspace(0.625*math.pi, 0.0, 20):
        x_vals.append(x_cen + 1.95*radius*math.cos(angle))
        y_vals.append(y_top + (0.2*height)*math.sin(angle))
    axis.plot(x_vals, y_vals, color=color)
    axis.plot([x_vals[-1]], [y_vals[-1]-arrow], 'v', 
                markerfacecolor=color, markeredgecolor='none', markersize=3)
    # Cross out the displaced cell
    #color = 'red'
    x_left, x_right = x1 + 6*radius, x1 + 8*radius
    y_bottom, y_top = y1 + 0.2*height + y2xscale*radius, y1 + 0.6*height - y2xscale*radius
    axis.plot([x_left, x_right], [y_top + y2xscale*radius, y_bottom - y2xscale*radius], color)
    axis.plot([x_left, x_right], [y_bottom - y2xscale*radius, y_top + y2xscale*radius], color)


def prettify_lineages(sub):
    radius = 0.25
    far_diff, near_diff = 0.15, 0.05
    y2xscale = 0.5/6
    text_left, text_diff = 1.5 - radius, far_diff + y2xscale*radius + 0.005
    cell = toolbox_schematic.Cell()
    cell.set_defaults(facecolor='#E6E6E6', edgecolor='0.5', y2xscale=y2xscale)
    cell.set_points([1.5, 1 + near_diff], [1.5, 1 + far_diff], radius)
    cell.draw(sub)
    sub.plot(cell.x_vals, cell.y1_vals, color='blue', lw=2)
    sub.text(text_left, 1 + text_diff, 'New pole lineage',
                                        color='blue', va='bottom', ha='left')
    cell.set_points([1.5, 1 - near_diff], [1.5, 1 - far_diff], radius)
    cell.draw(sub)
    sub.plot(cell.x_vals, cell.y2_vals, color='red', lw=2)
    sub.text(text_left, 1 - text_diff, 'Old pole lineage',
                                         color='red', va='top', ha='left')


def growth_finish(axis, toxic=False, xlim=0.25):
    ylim = 0.6
    draw_cell(axis, 0.1*xlim, 0.2*xlim, 0.1*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=toxic, arrow=0.01*xlim)
    draw_const_env(axis, 0.3*xlim, 0.00, 0.5*xlim, 0.2*ylim, y2xscale=ylim/xlim)
    # Finish off the axis
    axis.set_xlim([0.00, xlim])
    axis.set_ylim([0.00, ylim])


def glucose_finish(axis, toxic=False, xlim=0.25):
    ylim = 3.5
    draw_cell(axis, 0.1*xlim, 0.2*xlim, 0.1*ylim, 0.05*xlim, y2xscale=ylim/xlim, toxic=toxic, arrow=0.01*xlim)
    # Draw the chemostat
    draw_chemostat(axis, 0.3*xlim, 0.00, 0.5*xlim, 0.15*ylim)
    # Finish off the axis
    axis.plot([-1, 1],[3.24, 3.24], color='grey', linestyle='--')
    axis.text(0.0, 3.29, 'Inflow concentration', color='grey',
              ha = 'left', va='bottom')
    axis.set_xlim([0.00, xlim])
    axis.set_ylim([0.0, ylim])


def erjavec_finish(axis):
    axis.plot([-1, 4], [0, 0], color='grey', linestyle=':')
    axis.set_xlim([0, 3.2])
    axis.set_ylim([-0.05, 0.45])
