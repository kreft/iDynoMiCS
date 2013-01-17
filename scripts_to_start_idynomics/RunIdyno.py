#!/usr/bin/env python
# RuniDyno.py
# Author: Rob Clegg, Chinmay Kanchi
import sys, os, glob, time
from optparse import OptionParser
import __main__ as main


# Please set this to the directory where the iDynoMiCS java class files are
srcpath = 'set_to_source_directory'

if srcpath == 'set_to_source_directory':
    #get name of current file, so we can get the location in case -s is not provided
    myfilename = main.__file__
    mypath = os.path.dirname(os.path.abspath(myfilename))
    srcpath = os.path.join(os.path.split(mypath)[0],'bin')
    if (not os.path.exists(srcpath)):
        srcpath = os.path.join(os.path.split(mypath)[0],'src')


### Parse options
parser = OptionParser()

parser.add_option("-x", "--Xmin", dest="Xmin", default=300, 
        type="int", help="minimum memory allocation XMIN", metavar="XMIN")
parser.add_option("-X", "--Xmax", dest="Xmax", default=1000, 
        type="int",help="maximum memory allocation XMAX", metavar="XMAX")
parser.add_option("-m", "--multiples", dest="Multiples", default=1, 
        type="int", help="number of multiple runs of same protocol MULTIPLES", metavar="MULTIPLES")
parser.add_option("-s", "--src", dest="NPATH", default=srcpath, 
        help="source code directory", metavar="NPATH")

(options, args) = parser.parse_args()

### Set the source directory and check it exists
NPATH = os.path.abspath(options.NPATH)
if (not os.path.exists(NPATH)):
    print "Source directory %s doesn't exist!" %NPATH
    exit()
print "Using source at %s" %NPATH

### Set the classpath
CLASSPATH = ':'.join([jar for jar in glob.glob(os.path.join(NPATH, 'lib/*.jar'))])
CLASSPATH=NPATH+':'+CLASSPATH


### Compile the list of protocol files
PROTOCOLS = args
# we then check that this file exists
for protoFile in PROTOCOLS[:]:
    if protoFile.lower()[-4:] != '.xml':
        print 'All non-option arguments must be xml files, skipping %s'%(protoFile)
        PROTOCOLS.remove(protoFile)
        continue
    if (not os.path.exists(protoFile)):
        print "%s doesn't exist!" %protoFile
        PROTOCOLS.remove(protoFile)


### Run the simulations
if (len(PROTOCOLS)<1):
    print "No protocol file!"
else:
    for protoFile in PROTOCOLS:
        for m in range(options.Multiples):
            cmd = '''java -Xms%sm -Xmx%sm -cp %s idyno.Idynomics %s'''%(options.Xmin,options.Xmax,CLASSPATH,protoFile)
            os.system(cmd)
            time.sleep(60)
