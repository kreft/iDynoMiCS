#Stats script
import csv
import os
import toolbox_basic as basic
import toolbox_results as results
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy
from scipy import stats

def get_last_only(folder, f):
    os.chdir(folder+f+'/lastIter')
    output = results.AgentOutput(path='agent_Sum(last).xml')
    species = results.SpeciesOutput(output, 'OldieA')
    requirements={}
    cells = species.find_cells(requirements)
    t_a = (float(output.time))
    cell = cells[0]
    population_a = (float(cell.vars['population']))
    biomass_a = (float(cell.vars['mass']))
    growthrate_a = (float(cell.vars['growthRate']))
    species = results.SpeciesOutput(output, 'OldieB')
    requirements={}
    cells = species.find_cells(requirements)
    cell = cells[0]
    population_b = (float(cell.vars['population']))
    biomass_b = (float(cell.vars['mass']))
    growthrate_b = (float(cell.vars['growthRate']))
    if growthrate_a == 0:
        growthrate_b = abs(1/growthrate_b)
        gr = numpy.log(growthrate_b)
    elif growthrate_b == 0:
        growthrate_a = abs(growthrate_a/1)
        gr = numpy.log(growthrate_a)
    elif growthrate_a == 0 and growthrate_b == 0:
        gr = 0
    else:
        gr = abs(growthrate_a/growthrate_b)
        gr = numpy.log(gr)
    br = numpy.log(biomass_a/biomass_b)
    pr = numpy.log(population_a/population_b)
    ratio = [br, pr, gr]
    time = t_a
    return time, ratio
    
def get_all(folder, f):
    sim_dir_path = basic.check_path(folder+f+'/')
    file_dir = os.path.join(sim_dir_path, 'agent_Sum')
    basic.unzip_files(file_dir+'.zip')
    file_list = basic.file_list(file_dir)
    t_a, population_a, biomass_a, growthrate_a = [], [], [], []
    t_b, population_b, biomass_b, growthrate_b = [], [], [], []
    last_pop_a, last_biomass_a, last_growthrate_a = [], [], []
    last_pop_b, last_biomass_b, last_growthrate_b = [], [], []
    for filename in file_list:
        output = results.AgentOutput(path=filename)
        species = results.SpeciesOutput(output, 'OldieA')
        requirements={}
        cells = species.find_cells(requirements)
        #single_result = results.SingleResult()
        t_a.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population_a.append(float(cell.vars['population']))
        biomass_a.append(float(cell.vars['mass']))
        growthrate_a.append(float(cell.vars['growthRate']))
        species = results.SpeciesOutput(output, 'OldieB')
        requirements={}
        cells = species.find_cells(requirements)
        #single_result = results.SingleResult()
        t_b.append(float(output.time))
        if len(cells) == 0:
            continue
        cell = cells[0]
        population_b.append(float(cell.vars['population']))
        biomass_b.append(float(cell.vars['mass']))
        growthrate_b.append(float(cell.vars['growthRate']))
    basic.rm_dir(file_dir)
    lists_a, lists_b = [population_a, biomass_a, growthrate_a], [population_b, biomass_b, growthrate_b]
    t_a1, t_b1 = t_a, t_b
    for i in range(3):
        list1, list2, list3, list4 = t_a1, lists_a[i], t_b1, lists_b[i]
        t_a, lists_a[i] = (list(x) for x in zip(*sorted(zip(list1, list2), key=itemgetter(0))))
        t_b, lists_b[i] = (list(x) for x in zip(*sorted(zip(list3, list4), key=itemgetter(0))))
    biomass_ratio, population_ratio, growthrate_ratio = [], [], []
    for j in range(len(lists_a[1])):
        biomass_ratio.append(numpy.log(lists_a[1][j]/lists_b[1][j]))
        population_ratio.append(numpy.log(lists_a[0][j]/lists_b[0][j]))
        if lists_a[2][j] == 0 and lists_b[2][j] == 0:
            growthrate_ratio.append(0)
        elif lists_a[2][j] == 0:
            lists_b[2][j] = abs(1/(lists_b[2][j]))
            growthrate_ratio.append(numpy.log(lists_b[2][j]))
        elif lists_b[2][j] == 0:
            lists_a[2][j] = abs((lists_a[2][j])/1)
            growthrate_ratio.append(numpy.log(lists_a[2][j]))
        else:
            growthrate_ratio.append(abs(numpy.log(lists_a[2][j]/lists_b[2][j])))
        if j == len(lists_a[1])-1:
            last_pop_a.append(lists_a[0][j])
            last_pop_b.append(lists_b[0][j])
            last_biomass_a.append(lists_a[1][j])
            last_biomass_b.append(lists_b[1][j])
            last_growthrate_a.append(lists_a[2][j])
            last_growthrate_b.append(lists_b[2][j])
    return t_a, [biomass_ratio, population_ratio, growthrate_ratio]

def stats_test_folder(folder, end_only, save_stats, name):
    inside_folders = os.listdir(folder)
    fns, times, ratios = [], [], []
    inside_folders_new = []
    for a in inside_folders:
        if a[-4:] != '.csv' and a[-4:] != '.png' and a[-8:] != 'DS_Store':
            inside_folders_new.append(a)
    inside_folders = inside_folders_new
    for a in range(len(inside_folders)):
        new_files = os.listdir(folder+inside_folders[a])
        os.chdir(folder+inside_folders[a])
        for b in range(len(new_files)):
            c = b+1
            b = new_files[b]
            if b[-4:] != '.png' and b[-8:] != 'DS_Store':
                if end_only:
                    time, ratio = get_last_only(folder+inside_folders[a]+'/', b)
                    fns.append(inside_folders[a]+'_'+str(c))
                    times.append(time)
                    ratios.append(ratio)
                else:
                    time, ratio = get_all(folder+inside_folders[a]+'/', b)
                    fns.append(inside_folders[a]+'_'+str(c))
                    times.append(time)
                    ratios.append(ratio)
    if end_only:
        os.chdir(save_stats+'end_only/')
        with open(name+'_end.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Competition', 'End time', 'Biomass', 'Population size', 'Growthrate'])
            for a in range(len(fns)):
                writer.writerow([fns[a]]+[times[a]]+ratios[a])
    if not end_only:
        os.chdir(save_stats+'time_course/')
        with open(name+'_all.csv', 'w') as f:
            writer = csv.writer(f)
            for a in range(len(fns)):
                writer.writerow([fns[a]]+['Time']+times[a])
                writer.writerow([fns[a]]+['Biomass']+ratios[a][0])
                writer.writerow([fns[a]]+['Population']+ratios[a][1])
                writer.writerow([fns[a]]+['Growth rate']+ratios[a][2])
    return

def run_functions(save_stats, base_folder, end_only=True):
    folders = []
    for a in os.listdir(base_folder):
        if a != 'figures' and a[-4:] != '.csv' and a[-4:] != '.png' and a[-8:] != 'DS_Store':
            folders.append(a)
    
    for b in folders:
        folder, name = base_folder+b+'/', b
        print(folder)
        if folder[-4:] != '.png' and folder[-8:] != 'DS_Store':
            stats_test_folder(folder, end_only, save_stats, name)
    return