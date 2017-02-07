#!/usr/bin/python
#import numpy
import os
import sys
import toolbox_basic as basic
#import toolbox_results as results
#import toolbox_fitness as fitness
#import toolbox_idynomics
import random
import math

#this is for one species only

dir_path = basic.check_path(sys.argv[1])
file_path = os.path.join(dir_path, 'even_spacing_to_start.xml')
gridres = 4
density = 201
genealogy = 0.0
generation = 0.0
birthday = 0.0
activeBiomassGrowth = 300
activeBiomassRepair = 0.0
inactiveBiomassGrowth = 0.0
inactiveBiomassRepair = 0.0
growthRate = 0.0
volumeRate = 0.0
locationX = 0.0 #- or the same as the radius?
locationY = 0.000000 #+i = +4 = +4*(i-1)
locationZ = 0.000000
radius = 0.0
totalRadius = 0.0
birthX = 0.0
birthY = 0.0
birthZ = 0.0
age = 0.0
hasDied = 'false'
deathX = 0.0
deathY = 0.0
deathZ = 0.0
deathday = 0.0

script = ""

script += '<idynomics>\n'
script += '<simulation iterate="0" time="0.0" unit="hour">\n'
script += '<grid resolution="4" nI="65" nJ="65" nK="1"/>\n'
script += '<species name="OldieA" header="family,genealogy,generation,birthday,activeBiomassGrowth,activeBiomassRepair,inactiveBiomassGrowth,inactiveBiomassRepair,growthRate,volumeRate,locationX,locationY,locationZ,radius,totalRadius,birthX,birthY,birthZ,age,hasDied,deathX,deathY,deathZ,deathday">\n'
for i in range(16):
    locationY += gridres
    activeBiomassG = activeBiomassGrowth * math.pow(2, random.random())
    radius = math.sqrt(activeBiomassG/(density*gridres*math.pi))
    totalRadius = radius
    locationX = radius
    script += '%i,%i,%i,%i,%f,%f,%f,%f,%i,%i,%f,%f,%f,%f,%f,%i,%i,%i,%f,%s,%i,%i,%i,%i; \n'%(i,genealogy,generation,birthday,activeBiomassG,activeBiomassRepair,inactiveBiomassGrowth,inactiveBiomassRepair,growthRate,volumeRate,locationX,locationY,locationZ,radius,totalRadius,birthX,birthY,birthZ,age,hasDied,deathX,deathY,deathZ,deathday)
script += '</species> \n'
script += '</simulation> \n'
script += '</idynomics> \n'

with open(file_path, 'w') as f:
    f.write(script)
    
    
    
    
    
