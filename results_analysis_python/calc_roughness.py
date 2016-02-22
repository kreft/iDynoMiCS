#!/usr/bin/python
from __future__ import division
import numpy
import os
import sys
import toolbox_basic as basic
import toolbox_results as results
import calc_roughness_locations

results_folder = sys.argv[1]
starting_time = 48

# These should output a file ('cell_locations.xml') containing the locations of all cells at starting_time (columns 'X' and 'Y'; default 240 hours)
# as well as 'grid' = [nI, nJ, nK, res]
calc_roughness_locations.get_all_cells_location(results_folder)
grid = calc_roughness_locations.get_grid_size(results_folder)

# Need agents coordinates and domain size
# Matrix of agentgrid of size xdom/gridres x ydom/gridres - xdom and ydom are height and width of computational domain
# Convert agent coordinates to agentgrid elements - 1 signifies presence
# Define bulk and biofilm as separate bodies
# Identify biofilm front points as agent occupied elements with at least one bulk neighbour, set equivalent elements to 1 in equivalent matrix, frontgrid    
    

# Get domain info
nI = float(grid[0])
nJ = float(grid[1])
gridres = float(grid[3])

# create a grid of zeros
agentgrid = numpy.zeros((nI, nJ), dtype=int)

base_path = os.path.join('~', 'git', 'iDynoMiCS')
print base_path
base_path = basic.check_path(base_path)
print base_path
input_path = os.path.join(base_path, 'results', 'biofilms', 'Preliminary_3day', 'Roughness', sys.argv[1])
#output_path = os.path.join(input_path, 'figure_4_life_history_ABC.pdf')
paths = ['cell_locations.xml']
attributes = {'name':'locations', 'header':'X,Y'}

results_path = os.path.join(input_path, paths[0])
results_output = results.ResultsOutput(path=results_path)
result_set = results.ResultSet(results_output, attributes)
#set elements of the grid that have a cell in to have a value of 1 - all other grids should remain as 0
for r in result_set.members:
    x, y = float(r.vars['X']), float(r.vars['Y'])
    x, y = x/gridres, y/gridres
    x, y = int(x), int(y)
    if agentgrid[x, y] == 0:
            agentgrid[x, y] = 1  

#set front cells to have a value of 2, making sure that they aren't on the edges where the biofilm wraps around
for (i, j) in zip(range(nI), range(nJ)):
    if agentgrid[i][j] < 1:
            continue
    # Look at the level below, if we're not on the bottom
    if i > 0:
        if agentgrid[i-1][j] == 0:
            agentgrid[i][j] = 2
            continue
        # Always apply cyclic search for left/right in the level below
        if agentgrid[i-1][(j-1)%nJ] == 0:
            agentgrid[i][j] = 2
            continue
        if agentgrid[i-1][(j+1)%nJ] == 0:
            agentgrid[i][j] = 2
            continue
    # Look at the level above, if we're not at the top
    if i < nI - 1:
        if agentgrid[i+1][j] == 0:
            agentgrid[i][j] = 2
            continue
        # Always apply cyclic search for left/right in the level above
        if agentgrid[i+1][(j-1)%nJ] == 0:
            agentgrid[i][j] = 2
            continue
        if agentgrid[i+1][(j+1)%nJ] == 0:
            agentgrid[i][j] = 2
            continue
    # Always apply cyclic search for left/right in the same level
    if agentgrid[i][(j-1)%nJ] == 0:
        agentgrid[i][j] = 2
        continue
    if agentgrid[i][(j+1)%nJ] == 0:
        agentgrid[i][j] = 2
        continue

frontgrid = numpy.zeros((nI), dtype=int)
for (i, j) in zip(range(nI), range(nJ)):
    if agentgrid[i, j] == 2:
        frontgrid[i] += 1

Cfx = numpy.zeros((nI), dtype=int)
Pf = 0
num = 0
for i in range(0, nI):
    temp = frontgrid[i]/nJ
    Cfx[i] = temp
    Pf += temp
    num += i * temp
print ("Cfx = ", Cfx)
Xf = num/Pf

Sigmaf = 0
num2 = 0
for j in range(nI):
    num2 += abs(j-Xf)*Cfx[j]
Sigmaf = num2/Pf
Sigma = Sigmaf/Xf
print ("Sigmaf = ", Sigmaf, "File = ", sys.argv[1])
print ("Sigma = ", Sigma)
print ("Xf = ", Xf)
print ("Pf = ", Pf)