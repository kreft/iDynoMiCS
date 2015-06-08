#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import toolbox_basic

input_dir =  os.path.join("~", "Dropbox", "KreftGroup", "NRM")
input_dir = toolbox_basic.check_path(input_dir)

input_file = os.path.join(input_dir, "bac_000_070000.ml")
output_file = os.path.join(input_dir, "bac_000_070000.pov")

with open(input_file, 'Ur') as f:
    input_script = f.readlines()

species_names = {'1':'Eco', '2':'Ego'}

# TODO .incl files?
output_script = ""

for line in input_script:
    if line is '':
        continue
    columns = line.replace('\n', '').split('\t')
    output_script += 'sphere { <%s, %s, %s> %s texture { %s } }\n' \
        %(columns[2], columns[1], columns[3], columns[4], species_names[columns[0]])

with open(output_file, 'w') as f:
    f.write(output_script)

