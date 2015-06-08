#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import toolbox_basic

input_file = \
        os.path.join("~", "Dropbox", "KreftGroup", "NRM", "bac_000_070000.ml")

input_file = toolbox_basic.check_path(input_file)

with open(input_file, 'Ur') as f:
    input_script = f.readlines()

species_names = {1:'Eco', 2:'Ego'}

species_lists = [[], []]

input_column_heads = ['species', 'locationX', 'locationY', 'locationZ',
                      'radius', 'unknown1', 'genealogy', 'growthRate', 'unknown2']

for line in input_script:
    if line is '':
        continue
    columns = line.replace('\n', '').split('\t')
    '''
    for i in range(len(columns)):
        if '.' in columns[i]:
            columns[i] = float(columns[i])
        else:
            columns[i] = int(columns[i])
    '''
    species_id = int(columns[0])
    output_line = columns[1]
    
    

output_script = '''<idynomics>
 <simulation iterate="0" time="0.0" unit="hour">
<grid resolution="8.0" nI="33" nJ="33" nK="1"/>
'''

